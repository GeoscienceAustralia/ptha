! Compile with -DEVOLVE_TIMER to add detailed timing of the evolve loop
#ifdef EVOLVE_TIMER
#   define EVOLVE_TIMER_START(tname) call ds%evolve_timer%timer_start(tname)
#   define EVOLVE_TIMER_STOP(tname)  call ds%evolve_timer%timer_end(tname)
#else
#   define EVOLVE_TIMER_START(tname)
#   define EVOLVE_TIMER_STOP(tname)
#endif

! Use definition below to use jacobi iteration (needs more work to implement, but leaving unfinished code here for now)
!#define DISPERSIVE_JACOBI

! Use definition below to use spatially variable jacobi iteration. This skips
! iterations in areas with errors < threshold, with some overhead for the adaptivity.
!#define DISPERSIVE_JACOBI_ADAPTIVE

! Use definition below to use Peregrine dispersion which includes an additional topographic term
! (by default we ignore this term, similar to JAGURS and many other tsunami codes)
!#define DISPERSIVE_PEREGRINE

! FIXME: The tridiagonal solve method might have efficiency improved by precomputing
!     lower_uh, diag_uh, upper_uh, lower_vh, diag_vh, upper_vh
!     .... and related infor for RHS vh/uh terms....
! at the cost of more storage

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
    real(dp), parameter :: jacobi_overrelax = 0.0_dp !, tridiagonal_overrelax = 0.0_dp

    type dispersive_solver_type

        logical :: is_staggered_grid 
        real(dp), allocatable :: RHS(:,:,:), last_U(:,:,:)
            !! Right hand side, and value of U at the start of the (dispersive) timestep

#ifdef DISPERSIVE_JACOBI
        real(dp), allocatable :: offdiagonalAx(:,:,:), diagonalA(:,:,:)
            !! For jacobi iteration, work arrays used to separate diagonal and off-diagonal part of implicit terms

#else
        real(dp), allocatable :: Ax(:,:,:)
            !! For tridiagonal iteration, work array used to compute dispersive RHS terms
#endif

#ifdef REALFLOAT
        real(dp) :: tol = 1.0e-06_dp !! Solver tolerance (JACOBI ITERATION)
#else
        real(dp) :: tol = 1.0e-08_dp !! Solver tolerance (JACOBI ITERATION)
#endif
        integer(ip) :: last_iter = 0 !! Iteration count at convergence (JACOBI ITERATION)
        integer(ip) :: max_iter = 1000 !! Maximum number of iterations in a single solve call (JACOBI ITERATION)
        real(dp) :: max_err = -HUGE(1.0_dp) !! Error tracking (JACOBI ITERATION)

        real(dp) :: td1 = 0.0_dp, td2 = -1.0_dp 
            !! Depths (below MSL) used to linearly taper-off dispersion. 
            !! Set td1 > td2 to linearly reduce dispersive terms to zero between depths of td1 --> td2
            !! Default corresponds to no tapering. 

        integer(ip) :: tridiagonal_inner_iter = 2_ip
            !! The tridiagonal solver does this many iterations of 'x-sweep'+'ysweep' for each call to solve_tridiag

        type(quadratic_extrapolation_type) :: qet
            !! Allow quadratic extrapolation in time

        type(timer_type) :: evolve_timer

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
        !$OMP PARALLEL DO
        do j = 1, ny
            ds%RHS(:,j,UH:VH) = 0.0_dp
            ds%last_U(:,j,UH:ELV) = 0.0_dp
        end do
        !$OMP END PARALLEL DO

#ifdef DISPERSIVE_JACOBI
        write(log_output_unit, *) 'Solution with Jacobi iteration not implemented (most code below but needs a bit more work)'
        call generic_stop
        ! Workspace for jacobi iteration
        if(allocated(ds%offdiagonalAx)) deallocate(ds%offdiagonalAx)
        if(allocated(ds%diagonalA)) deallocate(ds%diagonalA)
        allocate(ds%offdiagonalAx(nx, ny, UH:VH), ds%diagonalA(nx, ny, UH:VH))
        !$OMP PARALLEL DO
        do j = 1, ny
            ds%offdiagonalAx(:,j,UH:VH) = 0.0_dp
            ds%diagonalA(:,j,UH:VH) = 1.0_dp
        end do
        !$OMP END PARALLEL DO
#else
        ! Workspace for tridiagonal iteration
        if(allocated(ds%Ax)) deallocate(ds%Ax)
        allocate(ds%Ax(nx,ny,UH:VH))
        !$OMP PARALLEL DO
        do j = 1, ny
            ds%Ax(:,j,UH:VH) = 0.0_dp
        end do
        !$OMP END PARALLEL DO
#endif

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

                solution_uh(1,j) = uh(1,j)
                solution_uh(size(uh,1),j) = uh(size(uh,1),j)
                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_iph_j*d_iph_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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

                solution_vh(1,j) = vh(1,j)
                solution_vh(size(vh,1), j) = vh(size(vh,1), j)
                do i = 2, size(vh, 1) - 1

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_jph*d_i_jph / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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

                solution_uh(1,j) = uh(1,j)
                solution_uh(size(uh,1),j) = uh(size(uh,1),j)
                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j)
                    d_i_j = merge(msl_linear - elev(i, j), 0.0_dp, elev(i,j) < msl_linear)

#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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

                solution_vh(1,j) = vh(1,j)
                solution_vh(size(vh,1), j) = vh(size(vh,1), j)
                do i = 2, size(vh, 1) - 1

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    !d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                    !    elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
                    d_i_j = merge(msl_linear - elev(i,j), 0.0_dp, elev(i,j) < msl_linear)

#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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

                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_iph_j*d_iph_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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
                !uh(:,j) = b_uh + (tridiagonal_overrelax) * (b_uh - uh(:,j))

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

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_jph*d_i_jph / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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
                !vh(i,:) = b_vh + (tridiagonal_overrelax) * (b_vh - vh(i,:))

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

                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    d_i_j = merge((msl_linear - elev(i,j)), 0.0_dp, elev(i,j) < msl_linear)

#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (R_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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
                !uh(:,j) = b_uh + (tridiagonal_overrelax) * (b_uh - uh(:,j))

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

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    d_i_j = merge((msl_linear - elev(i,j)), 0.0_dp, elev(i,j) < msl_linear)

#ifdef DISPERSIVE_PEREGRINE
                    !
                    ! Peregrine dispersion
                    !
                    dispersive_premult = d_i_j*d_i_j / (r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! Compute various inverse depth terms, which are zeroed at 'dry' cells.
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
                !vh(i,:) = b_vh + (tridiagonal_overrelax) * (b_vh - vh(i,:))

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
            ! Extrapolate forward in time
            call ds%qet%extrapolate_in_time(forward_time, U(:,:,UH:VH), &
                do_nothing_if_missing_times=.true., did_extrapolate=did_extrapolate)
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

        ds%last_iter = niter
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

!#ifdef DISPERSIVE_JACOBI
!    !
!    ! JACOBI ITERATION -- not yet integrated
!    !
!
!    subroutine linear_dispersive_staggered_matmult_JACOBI(elev, uh, vh, msl_linear, &
!            distance_bottom_edge, distance_left_edge, &
!            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, &
!            update_UH, update_VH)
!        !
!        ! Towards solving the linear dispersive equation as per:
!        !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
!        !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
!        ! Or
!        !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
!        !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
!        !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017,
!        !   111, 46–54
!        !
!        ! BACKGROUND:
!        !
!        ! The shallow water equation for depth-integrated-x-velocity, with a dispersive term, can be organised like:
!        !    d(uh)/dt + { other shallow water equations terms } = dispersive_term_uh
!        ! where dispersive_term_uh contains a time-derivative of some second-order spatial derivatives, see further explanation
!        ! in the comments at the start of this file.
!        !
!        ! This can be written as
!        !     uh  = [ shallow_water_solution ] + dt*dispersive_term_uh
!        ! Since dispersive_term_uh contains a single derivative wrt time, it turns out that (dt * dispersive_term_uh)
!        ! does not contain the term 'dt' (although it depends on the flow variables at the last timestep)
!        !
!        ! Suppose that ( dt * dispersive_term_uh) is written as:
!        !    (dt * dispersive_term_uh) = { A_uh%*%x - A_uh%*%x_lasttimestep }
!        ! , where 'A_uh' is a sparse matrix that discretizes the spatial derivatives in the dispersive term, 'x' is the flow
!        ! variables [stage, uh, vh, elevation] at the next time-step, and %*% is matrix multiplication.
!        !
!        ! Then this routine computes:
!        ! -- The diagonal part of A_uh for the uh variables, denoted D_uh (stored in diagonalA_uh).
!        !    Note this is NOT multiplied by x.
!        ! -- The off-diagonal part of A_uh when multiplied by x, i.e. [ (A_uh - D_uh)%*%x ] (stored in
!        !    offdiagonalAx_uh and offdiagonalAx_vh respectively).
!        !
!        ! All of the above is repeated for vh (with somewhat different equations in spherical coordinates).
!        !
!        ! This can be used to implement Jacobi iteration. The idea is that to solve:
!        !     uh = [shallow_water_solution_uh - A_uh%*%x_last_timestep ] + A_uh%*%x
!        !     vh = [shallow_water_solution_vh - A_vh%*%x_last_timestep ] + A_vh%*%x
!        ! we can iteratively solve (for i = 1, 2, 3, .....) :
!        !     uh_(i+1) = RHS_uh + ( A_uh - D_uh)%*%x_i + D_uh%*%x_(i+1)
!        !     vh_(i+1) = RHS_vh + ( A_vh - D_vh)%*%x_i + D_vh%*%x_(i+1)
!        ! until the differences between x_(i+1) and x_(i) are negligable. The RHS terms:
!        !     RHS_uh = [shallow_water_solution_uh - A_uh%*%x_last_timestep ]
!        !     RHS_vh = [shallow_water_solution_vh - A_vh%*%x_last_timestep ]
!        ! do not change in time.
!        !
!        !
!        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
!            dlon, dlat, msl_linear,  &
!            distance_bottom_edge(:), distance_left_edge(:)
!        real(dp), intent(inout) :: offdiagonalAx_uh(:,:), offdiagonalAx_vh(:,:), diagonalA_uh(:,:), diagonalA_vh(:,:)
!        logical, intent(in) :: update_UH, update_VH
!
!        integer(ip) :: i, j
!        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
!        real(dp) :: r_dlat, R_coslat_dlon_jp1, R_coslat_dlon_j, R_coslat_dlat_jp1, R_coslat_dlat_j, coslat_jp1h
!
!        ! BEWARE: THE JAGURS PAPER USES A DIFFERENT COORDINATE SYSTEM TO SWALS
!        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.
!
!        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
!        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, &
!        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)
!
!        if(update_UH) then
!            !$OMP DO
!            do j = 2, size(uh, 2) - 1
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
!                R_coslat_dlat   = R_coslat_dlon * dlat / dlon
!                coslat_jph = distance_bottom_edge(j+1) * dlat / (dlon * distance_left_edge(1))
!                coslat_jmh = distance_bottom_edge(j+0) * dlat / (dlon * distance_left_edge(1))
!
!                do i = 2, size(uh, 1) - 1
!                    ! UH dispersive term
!
!                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
!                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
!                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
!
!                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * R_coslat_dlon)
!
!                    ! The diagonal component of A_uh (related to uh(i,j) )
!                    diagonalA_uh(i,j) = dispersive_premult * 1.0_dp/R_coslat_dlon * (-2.0_dp)
!
!                    offdiagonalAx_uh(i,j) = &
!                        ! First subtract the diagonal_uh%*%x component
!                        -diagonalA_uh(i,j) * uh(i,j) + &
!                        ! Then include all terms in A_uh%*%x here (easier bookkeeping to write it this way)
!                        dispersive_premult * ( &
!                            ! 1/(R_coslat) * d/dlon(uh)
!                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
!                            ! 1/(R_coslat) * d/dlat( coslat * vh )
!                            1.0_dp / R_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
!                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) )
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        if(update_VH) then
!            !$OMP DO
!            do j = 2, size(vh, 2) - 1
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                r_dlat = distance_left_edge(1)
!                R_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
!                R_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
!                R_coslat_dlat_jp1 = R_coslat_dlon_jp1 * dlat / dlon
!                R_coslat_dlat_j   = R_coslat_dlon_j   * dlat / dlon
!
!                coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( r_dlat  * dlon )
!                coslat_jph  = distance_bottom_edge(j+1) * dlat / ( r_dlat  * dlon )
!                coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( r_dlat  * dlon )
!
!                do i = 2, size(vh, 1) - 1
!
!                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
!                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
!                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
!
!                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat)
!
!                    ! The diagonal component of A_vh (related to vh(i,j))
!                    diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp / R_coslat_dlat_jp1 * coslat_jph + &
!                                                               1.0_dp / R_coslat_dlat_j   * coslat_jph )
!
!                    offdiagonalAx_vh(i,j) = &
!                        ! First subtract the diagonal_vh%*%x component
!                        -diagonalA_vh(i,j) * vh(i,j) + &
!                        ! Then include all terms in A_vh%*%x here (easier bookkeeping to write it this way)
!                        dispersive_premult * ( &
!                            ! uh term
!                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
!                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) + &
!                            ! vh term
!                            ( 1.0_dp / R_coslat_dlat_jp1 * ( coslat_jp1h * vh(i, j+1) - coslat_jph * vh(i, j  )) - &
!                              1.0_dp / R_coslat_dlat_j   * ( coslat_jph  * vh(i, j  ) - coslat_jmh * vh(i, j-1)) ) )
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        !$OMP END PARALLEL
!
!    end subroutine
!
!    subroutine linear_dispersive_cellcentred_matmult_JACOBI(elev, uh, vh, msl_linear, &
!            distance_bottom_edge, distance_left_edge, &
!            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, &
!            update_UH, update_VH)
!        !
!        ! Towards solving the linear dispersive equation as per:
!        !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
!        !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
!        ! Or
!        !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
!        !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
!        !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017,
!        !   111, 46–54
!        !
!        ! BACKGROUND:
!        !   This is a cell-centred version of linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE
!        !   See that routine for documentation
!        !
!        !
!        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
!            dlon, dlat, msl_linear,  &
!            distance_bottom_edge(:), distance_left_edge(:)
!        real(dp), intent(inout) :: offdiagonalAx_uh(:,:), offdiagonalAx_vh(:,:), diagonalA_uh(:,:), diagonalA_vh(:,:)
!        logical, intent(in) :: update_UH, update_VH
!
!        integer(ip) :: i, j
!        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jp1, coslat_j, coslat_jm1, d_i_j, dispersive_premult
!        real(dp) :: r_dlat, R_coslat_dlon_jph, R_coslat_dlon_jmh, R_coslat_dlat_jph, R_coslat_dlat_jmh
!
!        ! BEWARE: THE JAGURS PAPER USES A DIFFERENT COORDINATE SYSTEM TO SWALS
!        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.
!
!        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
!        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, &
!        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)
!
!        if(update_UH) then
!            !!$OMP DO SCHEDULE(DYNAMIC, 10)
!            !!$OMP DO SCHEDULE(STATIC)
!            !$OMP DO
!            do j = 2, size(uh, 2) - 1
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  ))
!                R_coslat_dlat   = R_coslat_dlon * dlat / dlon
!
!                ! NOTE regarding array bounds. 
!                ! - We know j > 1 and j < size(uh,2). 
!                ! - Also distance_bottom_edge has size == size(uh,2)+1
!                coslat_jp1 = 0.5_dp * (distance_bottom_edge(j+2) + distance_bottom_edge(j+1)) * &
!                    dlat / (dlon * distance_left_edge(1))
!                coslat_jm1 = 0.5_dp * (distance_bottom_edge(j+0) + distance_bottom_edge(j-1)) * &
!                    dlat / (dlon * distance_left_edge(1))
!
!                do i = 2, size(uh, 1) - 1
!                    ! UH dispersive term
!
!                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
!                    !d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
!                    !    elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
!                    d_i_j = merge(msl_linear - elev(i, j), 0.0_dp, elev(i,j) < msl_linear)
!
!                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * R_coslat_dlon) !* &
!                        ! Linear taper
!                        !min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))
!
!                    ! The diagonal component of A_uh (related to uh(i,j) )
!                    diagonalA_uh(i,j) = dispersive_premult * 1.0_dp/R_coslat_dlon * (-2.0_dp)
!
!                    offdiagonalAx_uh(i,j) = &
!                        ! First subtract the diagonal_uh%*%x component
!                        -diagonalA_uh(i,j) * uh(i,j) + &
!                        ! Then include all terms in A_uh%*%x here (easier bookkeeping to write it this way)
!                        ! h0^2 / (3 R_coslat) * d/dlon
!                        dispersive_premult * ( &
!                            ! lon-diff of { 1/(R_coslat) * d/dlon(uh) }
!                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
!                            ! lon-diff of { 1/(R_coslat) * d/dlat( coslat * vh ) }
!                            1.0_dp / R_coslat_dlat * ( &
!                                (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
!                                (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp)  &
!                            )
!
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        if(update_VH) then
!            !!$OMP DO SCHEDULE(DYNAMIC, 10)
!            !!$OMP DO SCHEDULE(STATIC)
!            !$OMP DO
!            do j = 2, size(vh, 2) - 1
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                r_dlat = distance_left_edge(1)
!                R_coslat_dlon_jph = distance_bottom_edge(j+1)
!                R_coslat_dlon_jmh = distance_bottom_edge(j+0)
!                R_coslat_dlat_jph = R_coslat_dlon_jph * dlat / dlon
!                R_coslat_dlat_jmh = R_coslat_dlon_jmh * dlat / dlon
!
!                coslat_jp1 = 0.5_dp*(distance_bottom_edge(j+2)+distance_bottom_edge(j+1)) * dlat / (distance_left_edge(1) * dlon)
!                coslat_j   = 0.5_dp*(distance_bottom_edge(j+1)+distance_bottom_edge(j+0)) * dlat / (distance_left_edge(1) * dlon)
!                coslat_jm1 = 0.5_dp*(distance_bottom_edge(j+0)+distance_bottom_edge(j-1)) * dlat / (distance_left_edge(1) * dlon)
!
!
!                do i = 2, size(vh, 1) - 1
!
!                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
!                    !d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
!                    !    elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
!                    d_i_j = merge(msl_linear - elev(i,j), 0.0_dp, elev(i,j) < msl_linear)
!
!                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_dlat) !* &
!                        ! Linear taper
!                        !min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))
!
!                    ! The diagonal component of A_vh (related to vh(i,j))
!                    !diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp / R_coslat_dlat_jp1 * coslat_jph + &
!                    !                                           1.0_dp / R_coslat_dlat_j   * coslat_jph )
!                    diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp/R_coslat_dlat_jph * coslat_j  + &
!                                                               1.0_dp/R_coslat_dlat_jmh * coslat_j )
!
!                    offdiagonalAx_vh(i,j) = &
!                        ! First subtract the diagonal_vh%*%x component
!                        -diagonalA_vh(i,j) * vh(i,j) + &
!                        ! Then include all terms in A_vh%*%x here (easier bookkeeping to write it this way)
!                        ! h0**2 / (3 * R) * d/dlat * (
!                        dispersive_premult * ( &
!                            ! uh term
!                            ! 1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of central differences at j+1 and j.
!                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
!                            ! 1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of central differences at j and j-1
!                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) + &
!                            ! vh term
!                            ! 1/(R coslat ) d(coslat vh)/dlat @ i, j+1/2
!                            ( 1.0_dp / R_coslat_dlat_jph * ( coslat_jp1 * vh(i, j+1) - coslat_j   * vh(i, j  )) - &
!                            ! 1/(R coslat ) d(coslat vh)/dlat @ i, j-1/2
!                              1.0_dp / R_coslat_dlat_jmh * ( coslat_j   * vh(i, j  ) - coslat_jm1 * vh(i, j-1)) ) )
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        !$OMP END PARALLEL
!
!    end subroutine
!
!    subroutine linear_dispersive_solve_staggered_grid_JACOBI(ds, U, &
!            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear)
!        ! Use Jacobi iteration to solve the linear dispersive discretization presented in
!        ! Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
!        ! Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
!        class(dispersive_solver_type), intent(inout) :: ds
!        real(dp), intent(inout) :: U(:,:,:)
!            ! domain%U after an explicit shallow water update
!        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
!            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver
!
!        integer :: i, j, iter
!        real(dp) :: max_err, last_U
!
!        ! Get the explicit part of the dispersive terms
!        call linear_dispersive_staggered_matmult_JACOBI(&
!            ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
!            msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!            ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!            ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!            update_UH=.true., update_VH=.true.)
!
!        ! Setup the right-hand-side term, combining the shallow-water solution with the explicit part
!        ! of the dispersive term
!        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds, U)
!        !$OMP DO
!        do j = 2, size(U,2) - 1
!            ds%RHS(:,j,UH) = U(:,j,UH) - (ds%offdiagonalAx(:,j,UH) + ds%diagonalA(:,j,UH) * ds%last_U(:,j,UH))
!        end do
!        !$OMP END DO NOWAIT
!        ! Do-nothing boundary conditions
!        ds%RHS(:,1,UH) = U(:,1,UH)
!        ds%RHS(:,size(U,2),UH) = U(:,size(U,2),UH)
!        ds%RHS(1,:,UH) = U(1,:,UH)
!        ds%RHS(size(U,1),:,UH) = U(size(U,1),:,UH)
!
!        !$OMP DO
!        do j = 2, size(U,2) - 1
!            ds%RHS(:,j,VH) = U(:,j,VH) - (ds%offdiagonalAx(:,j,VH) + ds%diagonalA(:,j,VH) * ds%last_U(:,j,VH))
!        end do
!        !$OMP END DO
!        ! Do-nothing boundary conditions
!        ds%RHS(:,1,VH) = U(:,1,VH)
!        ds%RHS(:,size(U,2),VH) = U(:,size(U,2),VH)
!        ds%RHS(1,:,VH) = U(1,:,VH)
!        ds%RHS(size(U,1),:,VH) = U(size(U,1),:,VH)
!
!        !$OMP END PARALLEL
!
!        jacobi_iter: do iter = 1, ds%max_iter
!
!            ! Jacobi iteration
!            ds%last_iter = iter
!
!            !
!            !UH update
!            !
!            call linear_dispersive_staggered_matmult_JACOBI(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.true., update_VH=.false.)
!
!            max_err = 0.0_dp
!            !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
!            do j = 2, size(U,2) - 1
!                do i = 2, size(U, 1) - 1
!                    last_U = U(i,j,UH)
!                    U(i,j,UH) = (ds%RHS(i,j,UH) + ds%offdiagonalAx(i,j,UH))/(1.0_dp - ds%diagonalA(i,j,UH))
!                    U(i,j,UH) = U(i,j,UH) + jacobi_overrelax*(U(i,j,UH)-last_U)
!
!                    ! Record the max abs_uh_difference/depth, reducing to abs_uh_difference in depths < 1m
!                    max_err = max( max_err, &
!                        (abs(U(i,j,UH) - last_U)/&
!                         max(msl_linear - 0.5_dp * (U(i+1,j,ELV) + U(i,j,ELV)), 1.0_dp)) )
!                end do
!            end do
!            !$OMP END PARALLEL DO
!
!            !
!            ! VH update
!            !
!            call linear_dispersive_staggered_matmult_JACOBI(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.false., update_VH=.true.)
!
!            !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
!            do j = 2, size(U,2) - 1
!                do i = 2, size(U, 1) - 1
!                    last_U = U(i,j,VH)
!                    U(i,j,VH) = (ds%RHS(i,j,VH) + ds%offdiagonalAx(i,j,VH))/(1.0_dp - ds%diagonalA(i,j,VH))
!                    U(i,j,VH) = U(i,j,VH) + jacobi_overrelax*(U(i,j,VH)-last_U)
!                    ! Record the max abs_vh_difference/depth, reducing to abs_vh_difference in depths < 1m
!                    max_err = max( max_err, &
!                        (abs(U(i, j,VH) - last_U)/&
!                         max(msl_linear - 0.5_dp * (U(i,j,ELV) + U(i,j+1,ELV)), 1.0_dp)) )
!                end do
!            end do
!            !$OMP END PARALLEL DO
!
!            !print*, '      err', max_err, '; iter ', iter
!
!            ! Check for tolerance
!            if(max_err < ds%tol) exit jacobi_iter
!
!        end do jacobi_iter
!
!        !!print*, 'Jacobi iter: ', ds%last_iter, max_err!, maxval(abs(U(:,:,UH:VH) - ds%last_U(:,:,UH:VH)))
!        !if(ds%last_iter == ds%max_iter) then
!        !    write(log_output_unit, *) 'Jacobi iteration hit max iterations (', ds%max_iter, ') with error ', max_err
!        !end if
!
!    end subroutine
!
!    !subroutine linear_dispersive_solve_cellcentred_grid_JACOBI(ds, U, &
!    !        dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear)
!    !
!    !end subroutine
!
!#ifdef DISPERSIVE_JACOBI_ADAPTIVE
!
!    subroutine linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(elev, uh, vh, msl_linear, &
!            distance_bottom_edge, distance_left_edge, &
!            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, &
!            update_UH, update_VH, needs_update, i0, i1, td1, td2)
!        !
!        ! Towards solving the linear dispersive equation as per:
!        !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
!        !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
!        ! Or
!        !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
!        !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
!        !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017,
!        !   111, 46–54
!        !
!        ! BACKGROUND:
!        !
!        ! The shallow water equation for depth-integrated-x-velocity, with a dispersive term, can be organised like:
!        !    d(uh)/dt + { other shallow water equations terms } = dispersive_term_uh
!        ! where dispersive_term_uh contains a time-derivative of some second-order spatial derivatives, see further explanation
!        ! in the comments at the start of this file.
!        !
!        ! This can be written as
!        !     uh  = [ shallow_water_solution ] + dt*dispersive_term_uh
!        ! Since dispersive_term_uh contains a single derivative wrt time, it turns out that (dt * dispersive_term_uh)
!        ! does not contain the term 'dt' (although it depends on the flow variables at the last timestep)
!        !
!        ! Suppose that ( dt * dispersive_term_uh) is written as:
!        !    (dt * dispersive_term_uh) = { A_uh%*%x - A_uh%*%x_lasttimestep }
!        ! , where 'A_uh' is a sparse matrix that discretizes the spatial derivatives in the dispersive term, 'x' is the flow
!        ! variables [stage, uh, vh, elevation] at the next time-step, and %*% is matrix multiplication.
!        !
!        ! Then this routine computes:
!        ! -- The diagonal part of A_uh for the uh variables, denoted D_uh (stored in diagonalA_uh).
!        !    Note this is NOT multiplied by x.
!        ! -- The off-diagonal part of A_uh when multiplied by x, i.e. [ (A_uh - D_uh)%*%x ] (stored in
!        !    offdiagonalAx_uh and offdiagonalAx_vh respectively).
!        !
!        ! All of the above is repeated for vh (with somewhat different equations in spherical coordinates).
!        !
!        ! This can be used to implement Jacobi iteration. The idea is that to solve:
!        !     uh = [shallow_water_solution_uh - A_uh%*%x_last_timestep ] + A_uh%*%x
!        !     vh = [shallow_water_solution_vh - A_vh%*%x_last_timestep ] + A_vh%*%x
!        ! we can iteratively solve (for i = 1, 2, 3, .....) :
!        !     uh_(i+1) = RHS_uh + ( A_uh - D_uh)%*%x_i + D_uh%*%x_(i+1)
!        !     vh_(i+1) = RHS_vh + ( A_vh - D_vh)%*%x_i + D_vh%*%x_(i+1)
!        ! until the differences between x_(i+1) and x_(i) are negligable. The RHS terms:
!        !     RHS_uh = [shallow_water_solution_uh - A_uh%*%x_last_timestep ]
!        !     RHS_vh = [shallow_water_solution_vh - A_vh%*%x_last_timestep ]
!        ! do not change in time.
!        !
!        !
!        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
!            dlon, dlat, msl_linear,  &
!            distance_bottom_edge(:), distance_left_edge(:)
!        real(dp), intent(inout) :: offdiagonalAx_uh(:,:), offdiagonalAx_vh(:,:), diagonalA_uh(:,:), diagonalA_vh(:,:)
!        logical, intent(in) :: update_UH, update_VH, needs_update(:,:)
!        integer(ip), intent(in) :: i0(:), i1(:)
!        real(dp), intent(in) :: td1, td2 ! "Depths below msl" used to taper off dispersion
!
!        integer(ip) :: i, j
!        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
!        real(dp) :: r_dlat, R_coslat_dlon_jp1, R_coslat_dlon_j, R_coslat_dlat_jp1, R_coslat_dlat_j, coslat_jp1h
!
!        ! BEWARE: THE JAGURS PAPER USES A DIFFERENT COORDINATE SYSTEM TO SWALS
!        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.
!
!        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
!        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, needs_update, i0, i1, &
!        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH, td1, td2)
!
!        if(update_UH) then
!            !$OMP DO SCHEDULE(DYNAMIC, 10)
!            !!$OMP DO SCHEDULE(STATIC)
!            do j = 2, size(uh, 2) - 1
!                !if(.not. any(needs_update(:,j))) cycle
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
!                R_coslat_dlat   = R_coslat_dlon * dlat / dlon
!                coslat_jph = distance_bottom_edge(j+1) * dlat / (dlon * distance_left_edge(1))
!                coslat_jmh = distance_bottom_edge(j+0) * dlat / (dlon * distance_left_edge(1))
!
!                !do i = 2, size(uh, 1) - 1
!                do i = i0(j), i1(j)
!                    ! UH dispersive term
!                    if(.not. needs_update(i,j)) cycle
!
!                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
!                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
!                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
!
!                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * R_coslat_dlon) * &
!                        ! Linear taper
!                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))
!
!                    ! The diagonal component of A_uh (related to uh(i,j) )
!                    diagonalA_uh(i,j) = dispersive_premult * 1.0_dp/R_coslat_dlon * (-2.0_dp)
!
!                    offdiagonalAx_uh(i,j) = &
!                        ! First subtract the diagonal_uh%*%x component
!                        -diagonalA_uh(i,j) * uh(i,j) + &
!                        ! Then include all terms in A_uh%*%x here (easier bookkeeping to write it this way)
!                        dispersive_premult * ( &
!                            ! 1/(R_coslat) * d/dlon(uh)
!                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
!                            ! 1/(R_coslat) * d/dlat( coslat * vh )
!                            1.0_dp / R_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
!                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) )
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        if(update_VH) then
!            !$OMP DO SCHEDULE(DYNAMIC, 10)
!            !!$OMP DO SCHEDULE(STATIC)
!            do j = 2, size(vh, 2) - 1
!                !if(.not. any(needs_update(:,j))) cycle
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                r_dlat = distance_left_edge(1)
!                R_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
!                R_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
!                R_coslat_dlat_jp1 = R_coslat_dlon_jp1 * dlat / dlon
!                R_coslat_dlat_j   = R_coslat_dlon_j   * dlat / dlon
!
!                coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( r_dlat  * dlon )
!                coslat_jph  = distance_bottom_edge(j+1) * dlat / ( r_dlat  * dlon )
!                coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( r_dlat  * dlon )
!
!                !do i = 2, size(vh, 1) - 1
!                do i = i0(j), i1(j)
!                    if(.not. needs_update(i,j)) cycle
!
!                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
!                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
!                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
!
!                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat) * &
!                        ! Linear taper
!                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))
!
!
!                    ! The diagonal component of A_vh (related to vh(i,j))
!                    diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp / R_coslat_dlat_jp1 * coslat_jph + &
!                                                               1.0_dp / R_coslat_dlat_j   * coslat_jph )
!
!                    offdiagonalAx_vh(i,j) = &
!                        ! First subtract the diagonal_vh%*%x component
!                        -diagonalA_vh(i,j) * vh(i,j) + &
!                        ! Then include all terms in A_vh%*%x here (easier bookkeeping to write it this way)
!                        dispersive_premult * ( &
!                            ! uh term
!                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
!                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) + &
!                            ! vh term
!                            ( 1.0_dp / R_coslat_dlat_jp1 * ( coslat_jp1h * vh(i, j+1) - coslat_jph * vh(i, j  )) - &
!                              1.0_dp / R_coslat_dlat_j   * ( coslat_jph  * vh(i, j  ) - coslat_jmh * vh(i, j-1)) ) )
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        !$OMP END PARALLEL
!
!    end subroutine
!
!    subroutine linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(elev, uh, vh, msl_linear, &
!            distance_bottom_edge, distance_left_edge, &
!            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, &
!            update_UH, update_VH, needs_update, i0, i1, td1, td2)
!        !
!        ! Towards solving the linear dispersive equation as per:
!        !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
!        !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
!        ! Or
!        !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
!        !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
!        !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017,
!        !   111, 46–54
!        !
!        ! BACKGROUND:
!        !   This is a cell-centred version of linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE
!        !   See that routine for documentation
!        !
!        !
!        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
!            dlon, dlat, msl_linear,  &
!            distance_bottom_edge(:), distance_left_edge(:)
!        real(dp), intent(inout) :: offdiagonalAx_uh(:,:), offdiagonalAx_vh(:,:), diagonalA_uh(:,:), diagonalA_vh(:,:)
!        logical, intent(in) :: update_UH, update_VH, needs_update(:,:)
!        integer(ip), intent(in) :: i0(:), i1(:)
!        real(dp), intent(in) :: td1, td2 ! "Depth below MSL" used to linear taper off dispersive terms
!
!        integer(ip) :: i, j
!        real(dp) :: R_coslat_dlat, R_coslat_dlon, coslat_jp1, coslat_j, coslat_jm1, d_i_j, dispersive_premult
!        real(dp) :: r_dlat, R_coslat_dlon_jph, R_coslat_dlon_jmh, R_coslat_dlat_jph, R_coslat_dlat_jmh
!
!        ! BEWARE: THE JAGURS PAPER USES A DIFFERENT COORDINATE SYSTEM TO SWALS
!        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.
!
!        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
!        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, needs_update, i0, i1, &
!        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH, td1, td2)
!
!        if(update_UH) then
!            !$OMP DO SCHEDULE(DYNAMIC, 10)
!            !!$OMP DO SCHEDULE(STATIC)
!            do j = 2, size(uh, 2) - 1
!                !if(.not. any(needs_update(:,j))) cycle
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                R_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  ))
!                R_coslat_dlat   = R_coslat_dlon * dlat / dlon
!
!                ! NOTE regarding array bounds. 
!                ! - We know j > 1 and j < size(uh,2). 
!                ! - Also distance_bottom_edge has size == size(uh,2)+1
!                coslat_jp1 = 0.5_dp * (distance_bottom_edge(j+2) + distance_bottom_edge(j+1)) * &
!                    dlat / (dlon * distance_left_edge(1))
!                coslat_jm1 = 0.5_dp * (distance_bottom_edge(j+0) + distance_bottom_edge(j-1)) * &
!                    dlat / (dlon * distance_left_edge(1))
!
!                !do i = 2, size(uh, 1) - 1
!                do i = i0(j), i1(j)
!                    ! UH dispersive term
!                    if(.not. needs_update(i,j)) cycle
!
!                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
!                    !d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
!                    !    elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
!                    d_i_j = merge(msl_linear - elev(i, j), 0.0_dp, elev(i,j) < msl_linear)
!
!                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * R_coslat_dlon) * &
!                        ! Linear taper
!                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))
!
!                    ! The diagonal component of A_uh (related to uh(i,j) )
!                    diagonalA_uh(i,j) = dispersive_premult * 1.0_dp/R_coslat_dlon * (-2.0_dp)
!
!                    offdiagonalAx_uh(i,j) = &
!                        ! First subtract the diagonal_uh%*%x component
!                        -diagonalA_uh(i,j) * uh(i,j) + &
!                        ! Then include all terms in A_uh%*%x here (easier bookkeeping to write it this way)
!                        ! h0^2 / (3 R_coslat) * d/dlon
!                        dispersive_premult * ( &
!                            ! lon-diff of { 1/(R_coslat) * d/dlon(uh) }
!                            1.0_dp / R_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
!                            ! lon-diff of { 1/(R_coslat) * d/dlat( coslat * vh ) }
!                            1.0_dp / R_coslat_dlat * ( &
!                                (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
!                                (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp)  &
!                            )
!
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        if(update_VH) then
!            !$OMP DO SCHEDULE(DYNAMIC, 10)
!            !!$OMP DO SCHEDULE(STATIC)
!            do j = 2, size(vh, 2) - 1
!                !if(.not. any(needs_update(:,j))) cycle
!
!                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
!                ! cartesian coordinates.
!                r_dlat = distance_left_edge(1)
!                R_coslat_dlon_jph = distance_bottom_edge(j+1)
!                R_coslat_dlon_jmh = distance_bottom_edge(j+0)
!                R_coslat_dlat_jph = R_coslat_dlon_jph * dlat / dlon
!                R_coslat_dlat_jmh = R_coslat_dlon_jmh * dlat / dlon
!
!                coslat_jp1 = 0.5_dp*(distance_bottom_edge(j+2)+distance_bottom_edge(j+1)) * dlat / (distance_left_edge(1) * dlon)
!                coslat_j   = 0.5_dp*(distance_bottom_edge(j+1)+distance_bottom_edge(j+0)) * dlat / (distance_left_edge(1) * dlon)
!                coslat_jm1 = 0.5_dp*(distance_bottom_edge(j+0)+distance_bottom_edge(j-1)) * dlat / (distance_left_edge(1) * dlon)
!
!
!                !do i = 2, size(vh, 1) - 1
!                do i = i0(j), i1(j)
!                    if(.not. needs_update(i,j)) cycle
!
!                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
!                    !d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
!                    !    elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
!                    d_i_j = merge(msl_linear - elev(i,j), 0.0_dp, elev(i,j) < msl_linear)
!
!                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_dlat) * &
!                        ! Linear taper
!                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))
!
!                    ! The diagonal component of A_vh (related to vh(i,j))
!                    !diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp / R_coslat_dlat_jp1 * coslat_jph + &
!                    !                                           1.0_dp / R_coslat_dlat_j   * coslat_jph )
!                    diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp/R_coslat_dlat_jph * coslat_j  + &
!                                                               1.0_dp/R_coslat_dlat_jmh * coslat_j )
!
!                    offdiagonalAx_vh(i,j) = &
!                        ! First subtract the diagonal_vh%*%x component
!                        -diagonalA_vh(i,j) * vh(i,j) + &
!                        ! Then include all terms in A_vh%*%x here (easier bookkeeping to write it this way)
!                        ! h0**2 / (3 * R) * d/dlat * (
!                        dispersive_premult * ( &
!                            ! uh term
!                            ! 1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of central differences at j+1 and j.
!                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
!                            ! 1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of central differences at j and j-1
!                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) + &
!                            ! vh term
!                            ! 1/(R coslat ) d(coslat vh)/dlat @ i, j+1/2
!                            ( 1.0_dp / R_coslat_dlat_jph * ( coslat_jp1 * vh(i, j+1) - coslat_j   * vh(i, j  )) - &
!                            ! 1/(R coslat ) d(coslat vh)/dlat @ i, j-1/2
!                              1.0_dp / R_coslat_dlat_jmh * ( coslat_j   * vh(i, j  ) - coslat_jm1 * vh(i, j-1)) ) )
!                end do
!            end do
!            !$OMP END DO
!        end if
!
!        !$OMP END PARALLEL
!
!    end subroutine
!
!    subroutine linear_dispersive_solve_staggered_grid_JACOBI_ADAPTIVE(ds, U, &
!            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear, rhs_is_up_to_date, &
!            is_priority_domain_not_periodic, estimate_solution_forward_in_time, forward_time)
!        ! Use Jacobi iteration to solve the linear dispersive discretization presented in
!        ! Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
!        ! Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
!        !
!        ! This version is adaptive -- it skips computations where the error is (locally) small enough.
!        !
!        class(dispersive_solver_type), intent(inout) :: ds
!        real(dp), intent(inout) :: U(:,:,:)
!            ! domain%U after an explicit shallow water update
!        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
!            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver
!        logical, optional, intent(in) :: rhs_is_up_to_date
!        integer(ip), optional, intent(in) :: is_priority_domain_not_periodic(:,:) ! FIXME: Really optional?
!        logical, optional, intent(in) :: estimate_solution_forward_in_time
!        real(dp), optional, intent(in) :: forward_time
!
!
!        integer :: i, j, iter, i0, j0, i1, j1, ii0, ii1, l0, l1
!        real(dp) :: max_err, last_U, e1, e2
!        integer, parameter :: min_jacobi_iterations = 1, bw = 3, ob = 1
!        logical :: allowed_to_exit, rhs_ready, quadratic_extrap
!
!        rhs_ready = .FALSE.
!        if(present(rhs_is_up_to_date)) rhs_ready = rhs_is_up_to_date
!
!        if(present(estimate_solution_forward_in_time)) then
!            quadratic_extrap = .true.
!            if(.not. present(forward_time)) call generic_stop
!        end if
!
!        !$OMP PARALLEL DO DEFAULT(SHARED)
!        do j = 1, size(U,2)
!            ! Ensure that some iteration is performed
!            ds%needs_update(:,j) = .true.
!            ds%local_err(:,j) = 0.0_dp
!
!            ! i indices to update for each j -- initially look at all cells where we can compute central differences.
!            ds%i0(j) = 2
!            ds%i1(j) = size(U,1) - 1
!        end do
!        !$OMP END PARALLEL DO
!
!        if(.not. rhs_ready) then
!
!            ! Get the explicit part of the dispersive terms
!            call linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(&
!                ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.true., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
!                td1 = ds%td1, td2=ds%td2)
!
!            ! Setup the right-hand-side term, combining the shallow-water solution with the explicit part
!            ! of the dispersive term
!            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds, U)
!            !$OMP DO
!            do j = 2, size(U,2) - 1
!                ds%RHS(:,j,UH) = U(:,j,UH) - (ds%offdiagonalAx(:,j,UH) + ds%diagonalA(:,j,UH) * ds%last_U(:,j,UH))
!            end do
!            !$OMP END DO NOWAIT
!            !$OMP DO
!            do j = 2, size(U,2) - 1
!                ds%RHS(:,j,VH) = U(:,j,VH) - (ds%offdiagonalAx(:,j,VH) + ds%diagonalA(:,j,VH) * ds%last_U(:,j,VH))
!            end do
!            !$OMP END DO NOWAIT
!
!            !$OMP END PARALLEL
!        end if
!
!        ! FIXME: It would be nice to have the option to guess with a second-order-in-time extrapolation.
!        ! Then the solution will be second-order accurate at the first iteration, and maybe we don't need
!        ! as many iterations to make it good [although it will still depend on the wavelengths to be adjusted].
!        !call ds%store_last_U(U)
!        if(quadratic_extrap) call ds%qet%extrapolate_in_time(forward_time, U(:,:,UH:VH), do_nothing_if_missing_times=.TRUE.)
!
!        !allowed_to_exit = .FALSE. ! FIXME -- do we need this?
!
!        jacobi_iter: do iter = 1, ds%max_iter
!
!            ! Jacobi iteration
!            ds%last_iter = iter
!
!            !
!            !UH update
!            !
!            call linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.true., update_VH=.false., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
!                td1 = ds%td1, td2 = ds%td2)
!
!            max_err = 0.0_dp
!            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) &
!            !!$OMP PARALLEL DO SCHEDULE(STATIC) &
!            !$OMP DEFAULT(PRIVATE) SHARED(ds, U, msl_linear, is_priority_domain_not_periodic) &
!            !$OMP REDUCTION(max: max_err)
!            do j = 2, size(U,2) - 1
!                !do i = 2, size(U, 1) - 1
!                do i = ds%i0(j), ds%i1(j)
!                    if(.not. ds%needs_update(i,j)) cycle
!
!                    last_U = U(i,j,UH)
!                    U(i,j,UH) = (ds%RHS(i,j,UH) + ds%offdiagonalAx(i,j,UH))/(1.0_dp - ds%diagonalA(i,j,UH))
!                    U(i,j,UH) = U(i,j,UH) + jacobi_overrelax*(U(i,j,UH)-last_U)
!
!                    ! Velocity error metric
!                    e1 = abs(U(i,j,UH) - last_U)/&
!                         max(msl_linear - 0.5_dp * (U(i+1,j,ELV) + U(i,j,ELV)), 1.0_dp)
!                    ! Matrix-scale error metric
!                    !e2 = abs(U(i,j,UH)*(1.0_dp - ds%diagonalA(i,j,UH)) - (ds%RHS(i,j,UH) + ds%offdiagonalAX(i,j,UH)))
!
!                    ! Record the max abs_uh_difference/depth, reducing to abs_uh_difference in depths < 1m
!                    ds%local_err(i,j) = e1 !min(e1, e2) 
!                    if(present(is_priority_domain_not_periodic)) then
!                        if(is_priority_domain_not_periodic(i,j) > 0) max_err = max(max_err, ds%local_err(i,j))
!                    else
!                        max_err = max( max_err, ds%local_err(i,j) )
!                    end if
!                end do
!            end do
!            !$OMP END PARALLEL DO
!
!            !
!            ! VH update
!            !
!            call linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.false., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
!                td1 = ds%td1, td2 = ds%td2)
!
!            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) &
!            !!$OMP PARALLEL DO SCHEDULE(STATIC) &
!            !$OMP DEFAULT(PRIVATE) SHARED(ds, U, msl_linear, is_priority_domain_not_periodic) &
!            !$OMP REDUCTION(max: max_err)
!            do j = 2, size(U,2) - 1
!                !do i = 2, size(U, 1) - 1
!                do i = ds%i0(j), ds%i1(j)
!                    if(.not. ds%needs_update(i,j)) cycle
!
!                    last_U = U(i,j,VH)
!                    U(i,j,VH) = (ds%RHS(i,j,VH) + ds%offdiagonalAx(i,j,VH))/(1.0_dp - ds%diagonalA(i,j,VH))
!                    U(i,j,VH) = U(i,j,VH) + jacobi_overrelax*(U(i,j,VH)-last_U)
!
!                    ! Velocity error metric
!                    e1 = abs(U(i, j,VH) - last_U)/&
!                         max(msl_linear - 0.5_dp * (U(i,j,ELV) + U(i,j+1,ELV)), 1.0_dp)
!                    ! Matrix-scale error metric
!                    !e2 = abs(U(i,j,VH)*(1.0_dp - ds%diagonalA(i,j,VH)) - (ds%RHS(i,j,VH) + ds%offdiagonalAX(i,j,VH)))
!
!                    ! Record the max abs_vh_difference/depth, reducing to abs_vh_difference in depths < 1m
!                    ds%local_err(i,j) = max(ds%local_err(i,j), e1) !min(e1, e2))
!
!                    if(present(is_priority_domain_not_periodic)) then
!                        if(is_priority_domain_not_periodic(i,j) > 0) max_err = max(max_err, ds%local_err(i,j))
!                    else
!                        max_err = max( max_err, ds%local_err(i,j) )
!                    end if
!
!                end do
!            end do
!            !$OMP END PARALLEL DO
!
!            ds%max_err = max_err
!
!            ! Check for tolerance
!            if(ds%max_err < ds%tol) then
!                exit jacobi_iter
!                !if(allowed_to_exit) exit jacobi_iter
!                !! Once the error is less than the tolerance, we 
!                !! start doing full grid iterations again. 
!                !! Ideally we still meet the tolerance after one, and ensure 
!                !! the tolerance is met everywhere.
!                !! 
!                !allowed_to_exit = .true.
!                !ds%i0 = 2
!                !ds%i1 = size(U, 1)-1
!                !!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds)
!                !do j = 1, size(U,2)
!                !    ds%needs_update(:,j) = .true.
!                !end do
!                !!$OMP END PARALLEL DO
!            end if
!
!            !!if(iter > (min_jacobi_iterations - 1) .and. (.not. allowed_to_exit)) then
!            !    ! Determine which areas need updates.
!            !    ! We record cells if cells are 'near' a cell that has error > "some fraction of tol" (ds%needs_update(:,:)).
!            !    ! We also record, for each j, a range of i-indices that include all those needing updates.
!            !    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds)
!
!            !    !$OMP DO SCHEDULE(DYNAMIC, 10)
!            !    !!$OMP DO SCHEDULE(STATIC)
!            !    do j = 2, size(U,2) - 1
!            !        ! Update the 'needs_update' variable
!
!            !        ! We don't need to scan every i-index in the j'th row
!            !        ! Need to look look at range of active cells 'bw' rows above/below the current row,
!            !        ! as these might have error > tol. Need care to not go out of bounds. 
!            !        l0 = max(2, j-bw) ! safe j index below
!            !        l1 = min(size(U,2)-1, j+bw) ! safe j index above
!
!            !        ! i indices that we should look at in the j'th row
!            !        ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
!            !        ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above
!
!            !        do i = ii0 , ii1 
!            !            ! Check within local window for error > tol
!            !            i0 = max(1, i-bw)
!            !            i1 = min(size(U,1), i+bw)
!            !            j0 = max(1, j-bw)
!            !            j1 = min(size(U,2), j+bw)
!            !            ds%needs_update(i, j) = any(ds%local_err(i0:i1,j0:j1) > (ds%tol*ds%local_tol_factor))
!            !        end do
!            !    end do
!            !    !$OMP END DO
!
!                !!!$OMP DO SCHEDULE(DYNAMIC, 10)
!                !!$OMP DO SCHEDULE(STATIC)
!                !do j = 2, size(U, 2) - 1
!                !    ! Update the 'ds%i0, ds%i1' variables defining the i-range in which
!                !    ! we search for each j. Could not do this in the previous loop due to
!                !    ! race condition.
!
!                !    ! See the loop above for documentation of the next 4 lines
!                !    l0 = max(2, j-bw) ! safe j index below
!                !    l1 = min(size(U,2)-1, j+bw) ! safe j index above
!                !    ! FIXME: Race condition here? Considering update below
!                !    ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
!                !    ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above
!
!                !    ! To update in loop -- 
!                !    ! i0 = min index to update; initialise with a large value
!                !    ds%i0(j) = size(U, 1)
!                !    ! i1 = max index to update; initialise with a small value
!                !    ds%i1(j) = 1
!                !    do i = ii0, ii1
!                !        if(ds%needs_update(i,j)) then
!                !            ! Cell i,j is at least close to a cell that has error above tolerance.
!                !            ds%i0(j) = min(ds%i0(j), i)
!                !            ds%i1(j) = max(ds%i1(j), i)
!                !        end if
!                !    end do
!                !end do
!
!            !    !$OMP END PARALLEL
!            !end if
!
!        end do jacobi_iter
!
!        !!print*, 'Jacobi iter: ', ds%last_iter, max_err, '; shape(U) = ', shape(U) !count(ds%needs_update)*1.0_dp/(size(U,1)*size(U,2))
!        !if(ds%last_iter == ds%max_iter) then
!        !    write(log_output_unit, *) 'Jacobi iteration hit max iterations (', ds%max_iter, ') with error ', max_err
!        !end if
!
!    end subroutine
!
!    subroutine linear_dispersive_solve_cellcentred_JACOBI_ADAPTIVE(ds, U, &
!            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear, rhs_is_up_to_date, &
!            is_priority_domain_not_periodic, estimate_solution_forward_in_time, forward_time)
!        ! Use Jacobi iteration to solve the linear dispersive discretization presented in
!        ! Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
!        ! Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
!        !
!        ! Cellcentred version of linear_dispersive_solve_staggered_grid_JACOBI_ADAPTIVE
!        !
!        ! This version is adaptive -- it skips computations where the error is (locally) small enough.
!        !
!        class(dispersive_solver_type), intent(inout) :: ds
!        real(dp), intent(inout) :: U(:,:,:)
!            ! domain%U after an explicit shallow water update
!        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
!            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver
!        logical, optional, intent(in) :: rhs_is_up_to_date
!        integer(ip), optional, intent(in) :: is_priority_domain_not_periodic(:,:) ! FIXME: Currently not optional, but should it be?!
!        logical, optional, intent(in) :: estimate_solution_forward_in_time
!        real(dp), optional, intent(in) :: forward_time
!
!        integer(ip) :: i, j, iter, i0, j0, i1, j1, ii0, ii1, l0, l1
!        real(dp) :: max_err, last_U
!        integer(ip), parameter :: bw = 3
!        logical :: rhs_ready, quadratic_extrap
!
!        rhs_ready = .FALSE.
!        if(present(rhs_is_up_to_date)) rhs_ready = rhs_is_up_to_date
!
!        quadratic_extrap = .FALSE.
!        if(present(estimate_solution_forward_in_time)) then
!            quadratic_extrap = .true.
!            if(.not. present(forward_time)) call generic_stop
!        end if
!
!        !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i, i0, i1, j0, j1)
!        do j = 1, size(U,2)
!            ! Ensure that some iteration is performed
!            !ds%needs_update(:,j) = .true.
!            ds%needs_update(:,j) = (is_priority_domain_not_periodic(:,j) > 0_ip)
!
!            !-! Try buffering the 'needs update'
!            !do i = 1, size(U,1)
!            !    i0 = max(i-1,1)
!            !    i1 = min(i+1,size(U,1))
!            !    j0 = max(j-1,1)
!            !    j1 = min(j+1,size(U,2))
!            !    ds%needs_update(i,j) = any(is_priority_domain_not_periodic(i0:i1,j0:j1) > 0_ip) !.true.
!            !end do
!
!            ds%local_err(:,j) = 0.0_dp
!            ! i indices to update for each j -- initially look at all cells where we can compute central differences.
!            ds%i0(j) = 2
!            ds%i1(j) = size(U,1) - 1
!        end do
!        !$OMP END PARALLEL DO
!
!        if(.not. rhs_ready) then
!
!            ! Get the explicit part of the dispersive terms
!            call linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(&
!                ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.true., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
!                td1 = ds%td1, td2 = ds%td2)
!
!            ! Setup the right-hand-side term, combining the shallow-water solution with the explicit part
!            ! of the dispersive term
!            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds, U)
!            !$OMP DO
!            do j = 2, size(U,2) - 1
!                ds%RHS(:,j,UH) = U(:,j,UH) - (ds%offdiagonalAx(:,j,UH) + ds%diagonalA(:,j,UH) * ds%last_U(:,j,UH))
!            end do
!            !$OMP END DO NOWAIT
!            !$OMP DO
!            do j = 2, size(U,2) - 1
!                ds%RHS(:,j,VH) = U(:,j,VH) - (ds%offdiagonalAx(:,j,VH) + ds%diagonalA(:,j,VH) * ds%last_U(:,j,VH))
!            end do
!            !$OMP END DO NOWAIT
!
!            !$OMP END PARALLEL
!        end if
!
!        ! FIXME: It would be nice to have the option to guess with a second-order-in-time extrapolation.
!        ! Then the solution will be second-order accurate at the first iteration, and maybe we don't need
!        ! many iterations to make it good [although it will still depend on the wavelengths to be adjusted].
!        !call ds%store_last_U(U)
!        if(quadratic_extrap) call ds%qet%extrapolate_in_time(forward_time, U(:,:,UH:VH), do_nothing_if_missing_times=.TRUE.)
!
!        jacobi_iter: do iter = 1, ds%max_iter
!
!            ! Jacobi iteration
!            ds%last_iter = iter
!
!            !
!            !UH update
!            !
!            call linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.true., update_VH=.false., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
!                td1 = ds%td1, td2 = ds%td2)
!
!            max_err = 0.0_dp
!            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) &
!            !!$OMP PARALLEL DO SCHEDULE(STATIC) &
!            !$OMP DEFAULT(PRIVATE) SHARED(ds, U, msl_linear, is_priority_domain_not_periodic) &
!            !$OMP REDUCTION(max: max_err)
!            do j = 2, size(U,2) - 1
!                !do i = 2, size(U, 1) - 1
!                do i = ds%i0(j), ds%i1(j)
!                    if(.not. ds%needs_update(i,j)) cycle
!
!                    last_U = U(i,j,UH)
!                    U(i,j,UH) = (ds%RHS(i,j,UH) + ds%offdiagonalAx(i,j,UH))/(1.0_dp - ds%diagonalA(i,j,UH))
!                    U(i,j,UH) = U(i,j,UH) + jacobi_overrelax*(U(i,j,UH)-last_U)
!
!                    ! Record the max abs_uh_difference/depth, reducing to abs_uh_difference in depths < 1m
!                    ds%local_err(i,j) = ( abs(U(i,j,UH) - last_U)/&
!                         max(msl_linear - U(i,j,ELV), 1.0_dp) )
!
!                    if(present(is_priority_domain_not_periodic)) then
!                        if(is_priority_domain_not_periodic(i,j) > 0) max_err = max(max_err, ds%local_err(i,j))
!                    else
!                        max_err = max( max_err, ds%local_err(i,j) )
!                    end if
!                end do
!            end do
!            !$OMP END PARALLEL DO
!
!            !
!            ! VH update
!            !
!            call linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
!                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
!                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
!                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
!                update_UH=.false., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
!                td1 = ds%td1, td2 = ds%td2)
!
!            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) &
!            !!$OMP PARALLEL DO SCHEDULE(STATIC) &
!            !$OMP DEFAULT(PRIVATE) SHARED(ds, U, msl_linear, is_priority_domain_not_periodic) &
!            !$OMP REDUCTION(max: max_err)
!            do j = 2, size(U,2) - 1
!                !do i = 2, size(U, 1) - 1
!                do i = ds%i0(j), ds%i1(j)
!                    if(.not. ds%needs_update(i,j)) cycle
!
!                    last_U = U(i,j,VH)
!                    U(i,j,VH) = (ds%RHS(i,j,VH) + ds%offdiagonalAx(i,j,VH))/(1.0_dp - ds%diagonalA(i,j,VH))
!                    U(i,j,VH) = U(i,j,VH) + jacobi_overrelax*(U(i,j,VH)-last_U)
!
!                    ! Record the max abs_vh_difference/depth, reducing to abs_vh_difference in depths < 1m
!                    ds%local_err(i,j) = max(ds%local_err(i,j), abs(U(i, j,VH) - last_U)/&
!                         max(msl_linear - U(i,j,ELV), 1.0_dp))
!
!                    if(present(is_priority_domain_not_periodic)) then
!                        if(is_priority_domain_not_periodic(i,j) > 0) max_err = max(max_err, ds%local_err(i,j))
!                    else
!                        max_err = max( max_err, ds%local_err(i,j) )
!                    end if
!                end do
!            end do
!            !$OMP END PARALLEL DO
!
!            ds%max_err = max_err
!
!            ! Check for tolerance
!            if(ds%max_err < ds%tol) then
!                exit jacobi_iter
!            end if
!
!            !if(iter > (min_jacobi_iterations - 1) .and. (.not. allowed_to_exit)) then
!                !! Determine which areas need updates.
!                !! We record cells if cells are 'near' a cell that has error > "some fraction of tol" (ds%needs_update(:,:)).
!                !! We also record, for each j, a range of i-indices that include all those needing updates.
!                !!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds)
!
!                !!!$OMP DO SCHEDULE(DYNAMIC, 10)
!                !!$OMP DO SCHEDULE(STATIC)
!                !do j = 2, size(U,2) - 1
!                !    ! Update the 'needs_update' variable
!
!                !    ! We don't need to scan every i-index in the j'th row
!                !    ! Need to look look at range of active cells 'bw' rows above/below the current row,
!                !    ! as these might have error > tol. Need care to not go out of bounds. 
!                !    l0 = max(2, j-bw) ! safe j index below
!                !    l1 = min(size(U,2)-1, j+bw) ! safe j index above
!
!                !    ! i indices that we should look at in the j'th row
!                !    ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
!                !    ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above
!
!                !    do i = ii0 , ii1 
!                !        ! Check within local window for error > tol
!                !        i0 = max(1, i-bw)
!                !        i1 = min(size(U,1), i+bw)
!                !        j0 = max(1, j-bw)
!                !        j1 = min(size(U,2), j+bw)
!                !        ds%needs_update(i, j) = any(ds%local_err(i0:i1,j0:j1) > (ds%tol*ds%local_tol_factor))
!                !    end do
!                !end do
!                !!$OMP END DO
!
!                !!!$OMP DO SCHEDULE(DYNAMIC, 10)
!                !!$OMP DO SCHEDULE(STATIC)
!                !do j = 2, size(U, 2) - 1
!                !    ! Update the 'ds%i0, ds%i1' variables defining the i-range in which
!                !    ! we search for each j. Could not do this in the previous loop due to
!                !    ! race condition.
!
!                !    ! See the loop above for documentation of the next 4 lines
!                !    l0 = max(2, j-bw) ! safe j index below
!                !    l1 = min(size(U,2)-1, j+bw) ! safe j index above
!                !    ! FIXME: Race condition here -- see update below
!                !    ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
!                !    ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above
!
!                !    ! To update in loop -- 
!                !    ! i0 = min index to update; initialise with a large value
!                !    ds%i0(j) = size(U, 1)
!                !    ! i1 = max index to update; initialise with a small value
!                !    ds%i1(j) = 1
!                !    do i = ii0, ii1
!                !        if(ds%needs_update(i,j)) then
!                !            ! Cell i,j is at least close to a cell that has error above tolerance.
!                !            ds%i0(j) = min(ds%i0(j), i)
!                !            ds%i1(j) = max(ds%i1(j), i)
!                !        end if
!                !    end do
!                !end do
!                !!$OMP END DO
!
!                !!$OMP END PARALLEL
!            !end if
!
!        end do jacobi_iter
!
!        !print*, 'Jacobi iter: ', ds%last_iter, max_err, ', shape(U): ', shape(U) !, count(ds%needs_update)*1.0_dp/(size(U,1)*size(U,2))
!        !if(ds%last_iter == ds%max_iter) then
!        !    write(log_output_unit, *) 'Jacobi iteration hit max iterations (', ds%max_iter, ') with error ', max_err
!        !end if
!
!    end subroutine
!#endif
!#endif
end module
