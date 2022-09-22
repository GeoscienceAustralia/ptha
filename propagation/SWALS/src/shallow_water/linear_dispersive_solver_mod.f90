! Multigrid ideas might help
! https://petsc.org/release/src/ksp/ksp/tutorials/ex22f.F90.html

! Use definition below to use jacobi iteration
#define DISPERSIVE_JACOBI

! Use definition below to use iterative tridiagonal solve
!#define DISPERSIVE_TRIDIAG

module linear_dispersive_solver_mod
    !
    ! Type for solving the linear dispersive equation as per:
    !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
    !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
    ! Or
    !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
    !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
    !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017, 111,
    !   46–54
    ! The equations are:
    !     d( uh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !        h0^2 / ( 3 R_coslat) * d/dlon [ 1/R_coslat ( d^2( uh ) / (dt dlon) + d^2(vh coslat) / (dt dlat) ) ]
    !     d (vh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !        h0^2 / ( 3 R       ) * d/dlat [ 1/R_coslat ( d^2( uh ) / (dt dlon) + d^2(vh coslat) / (dt dlat) ) ]
    !
    ! Where 'h0' is the depth below MSL_linear. Note that Baba et al use 'co-latitude' instead of latitude, so my
    ! cos(theta) is their sin(theta)
    !
    ! ** Notice we can bring the time-derivative out to the front. This is key to solving the equation numerically. **
    !     d( uh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !        d/dt * { h0^2 / ( 3 R_coslat) * d/dlon [ 1/R_coslat ( d( uh ) / (dlon) + d(vh coslat) / (dlat) ) ] }
    !     d (vh ) / dt + {...nonlinear-shallow-water-terms...} = &
    !        d/dt * { h0^2 / ( 3 R       ) * d/dlat [ 1/R_coslat ( d( uh ) / (dlon) + d(vh coslat) / (dlat) ) ] }
    !
    ! To implement this in both cartesian and spherical coordinates, note that
    !     radius_earth = distance_cell_left_edge / dy_in_radians ! This is independent of the latitude or longitude.
    ! and
    !     R_coslat = distance_cell_bottom_edge / dx_in_radians ! Need to compute the cell-bottom distance at the right 'latitude'.
    ! So
    !     coslat = (distance_cell_bottom_edge / dx_in_radians) / ( distance_cell_left_edge / dy_in_radians )
    !     ! Note the factor to convert from degrees to radians will drop out in this equation.
    ! Also
    !     R_coslat_dlon = distance_cell_bottom_edge ! Here we need to compute the cell-bottom distance at the right 'latitude'.
    ! and
    !     R_dlat = distance_cell_left_edge ! This is independent of the latitude or longitude
    ! and
    !     R_dlat_coslat = distance_cell_left_edge * coslat = ( distance_cell_bottom_edge / dx_in_radians) * dy_in_radians


    !
    ! type(dispersive_solver_type) :: ds
    ! ! Before timestepping
    ! call ds%setup(nx, ny) ! [nx,ny] = domain%nx
    !
    ! timestepping_loop: do
    !     ...
    !     ! At start of timestep, we must backup domain%U, because the
    !     ! dispersive solver will need to know this value (as it has a time derivative).
    !     call ds%store_last_U(domain%U)
    !     ...
    !     ! Evolve one shallow water timestep
    !     ...
    !     ! Add in the dispersive term
    !     call ds%solve_staggered_grid(&
    !         domain%U, domain%dx(1), domain%dx(2), &
    !         domain%distance_bottom_edge, domain%distance_left_edge, &
    !         domain%msl_linear)
    !     ! Check for convergence
    !     if(ds%last_iter == ds%max_iter) then
    !         ! note convergence problems here
    !     end if
    !     ...
    ! end do timestepping_loop

    use global_mod, only: dp, ip
    use logging_mod, only: log_output_unit

    implicit none

    integer, parameter :: STG=1, UH=2, VH=3, ELV=4
    real(dp), parameter :: jacobi_overrelax = 0.0_dp, tridiagonal_overrelax = 0.3_dp

    type dispersive_solver_type
        !! Holds work arrays
        real(dp), allocatable :: RHS(:,:,:), last_U(:,:,:)
        real(dp), allocatable :: offdiagonalAx(:,:,:), diagonalA(:,:,:)

#ifdef REALFLOAT
        real(dp) :: tol = 1.0e-06_dp ! Solver tolerance -- iteration until this accuracy is reached
#else
        real(dp) :: tol = 1.0e-08_dp ! Solver tolerance -- iteration until this accuracy is reached
#endif
        integer(ip) :: max_iter = 2000 ! Maximum number of iterations
        integer(ip) :: last_iter = 0 ! Iteration count at convergence
        contains
        procedure :: store_last_U => store_last_U
        procedure :: setup => setup_dispersive_solver

#ifdef DISPERSIVE_JACOBI
        procedure :: solve_staggered_grid => dispersive_solve_staggered_grid_JACOBI
#endif
#ifdef DISPERSIVE_TRIDIAG
        procedure :: solve_staggered_grid => dispersive_solve_staggered_grid_TRIDIAG
#endif

    end type

    contains

    subroutine store_last_U(ds, U)
        ! For timestepping we need to retain U prior to the shallow water timestep.
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(in) :: U(:,:,:)

        integer(ip) :: j
        !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds, U)
        do j = 1, size(U,2)
            ds%last_U(:,j,UH:ELV) = U(:,j,UH:ELV)
        end do
        !$OMP END PARALLEL DO

    end subroutine

    subroutine setup_dispersive_solver(ds, nx, ny)
        !! Allocate workspace
        class(dispersive_solver_type), intent(inout) :: ds
        integer(ip), intent(in) :: nx, ny

        integer(ip) :: i, j

        if(allocated(ds%RHS)) deallocate(ds%RHS)
        if(allocated(ds%last_U)) deallocate(ds%last_U)

        if(allocated(ds%offdiagonalAx)) deallocate(ds%offdiagonalAx)
        if(allocated(ds%diagonalA)) deallocate(ds%diagonalA)

        allocate(ds%RHS(nx, ny, UH:VH), ds%last_U(nx, ny, UH:ELV) &
            , ds%offdiagonalAx(nx, ny, UH:VH), ds%diagonalA(nx, ny, UH:VH))

        ! OMP first-touch allocation
        !$OMP PARALLEL DO
        do j = 1, ny
            ds%RHS(:,j,UH:VH) = 0.0_dp
            ds%last_U(:,j,UH:ELV) = 0.0_dp

            ds%offdiagonalAx(:,j,UH:VH) = 0.0_dp
            ds%diagonalA(:,j,UH:VH) = 0.0_dp
        end do
        !$OMP END PARALLEL DO

    end subroutine


    !
    !
    !
    subroutine linear_dispersive_matmult(elev, uh, vh, msl_linear, distance_bottom_edge, distance_left_edge, &
            dlon, dlat, solution_uh, solution_vh, update_UH, update_VH)
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
        ! Then this routine computes (uh - A_uh%*%x) -- note the leading 'uh' accounts for the LHS uh term, and we have moved
        ! A_uh%*%x to the LHS.
        !
        ! The above is repeated for vh (with somewhat different equations in spherical coordinates).
        !
        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
            dlon, dlat, msl_linear,  &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(out) :: solution_uh(:,:), solution_vh(:,:)
        logical, intent(in) :: update_UH, update_VH

        integer(ip) :: i, j
        real(dp) :: r_coslat_dlat, r_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
        real(dp) :: r_dlat, r_coslat_dlon_jp1, r_coslat_dlon_j, r_coslat_dlat_jp1, r_coslat_dlat_j, coslat_jp1h

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, solution_uh, solution_vh, &
        !$OMP                                  dlat, dlon, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 1, size(uh, 2)

                if(j == 1 .or. j == size(uh, 2)) then
                    solution_uh(:,j) = uh(:,j)
                    cycle
                end if

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
                r_coslat_dlat   = r_coslat_dlon * dlat / dlon
                coslat_jph = distance_bottom_edge(j+1) * dlat / (dlon * distance_left_edge(1))
                coslat_jmh = distance_bottom_edge(j+0) * dlat / (dlon * distance_left_edge(1))

                solution_uh(1,j) = uh(1,j)
                solution_uh(size(uh,1),j) = uh(size(uh,1),j)
                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)

                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * r_coslat_dlon)

                    solution_uh(i,j) = uh(i,j) - &
                        ! Include all terms in A_uh%*%x here
                        dispersive_premult * ( &
                            ! 1/(R_coslat) * d/dlon(uh)
                            1.0_dp / r_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
                            ! 1/(R_coslat) * d/dlat( coslat * vh )
                            1.0_dp / r_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) )
                end do
            end do
            !$OMP END DO NOWAIT
        end if

        if(update_VH) then
            !$OMP DO
            do j = 2, size(vh, 2) - 1

                if(j == 1 .or. j == size(vh, 2)) then
                    solution_vh(:,j) = vh(:,j)
                    cycle
                end if

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_dlat = distance_left_edge(1)
                r_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
                r_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
                r_coslat_dlat_jp1 = r_coslat_dlon_jp1 * dlat / dlon
                r_coslat_dlat_j   = r_coslat_dlon_j   * dlat / dlon

                coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( distance_left_edge(1)  * dlon )
                coslat_jph  = distance_bottom_edge(j+1) * dlat / ( distance_left_edge(1)  * dlon )
                coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( distance_left_edge(1)  * dlon )

                solution_vh(1,j) = vh(1,j)
                solution_vh(size(vh,1), j) = vh(size(vh,1), j)
                do i = 2, size(vh, 1) - 1

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)

                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat)

                    solution_vh(i,j) = vh(i,j) - &
                        ! Include all terms in A_vh%*%x here
                        dispersive_premult * ( &
                            ! uh term
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) + &
                            ! vh term
                            ( 1.0_dp / r_coslat_dlat_jp1 * ( coslat_jp1h * vh(i, j+1) - coslat_jph * vh(i, j  )) - &
                              1.0_dp / r_coslat_dlat_j   * ( coslat_jph  * vh(i, j  ) - coslat_jmh * vh(i, j-1)) ) )
                end do
            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine

    !
    ! JACOBI ITERATION
    !

    subroutine linear_dispersive_matmult_JACOBI(elev, uh, vh, msl_linear, distance_bottom_edge, distance_left_edge, &
            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, update_UH, update_VH)
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
        ! Then this routine computes:
        ! -- The diagonal part of A_uh for the uh variables, denoted D_uh (stored in diagonalA_uh).
        !    Note this is NOT multiplied by x.
        ! -- The off-diagonal part of A_uh when multiplied by x, i.e. [ (A_uh - D_uh)%*%x ] (stored in
        !    offdiagonalAx_uh and offdiagonalAx_vh respectively).
        !
        ! All of the above is repeated for vh (with somewhat different equations in spherical coordinates).
        !
        ! This can be used to implement Jacobi iteration. The idea is that to solve:
        !     uh = [shallow_water_solution_uh - A_uh%*%x_last_timestep ] + A_uh%*%x
        !     vh = [shallow_water_solution_vh - A_vh%*%x_last_timestep ] + A_vh%*%x
        ! we can iteratively solve (for i = 1, 2, 3, .....) :
        !     uh_(i+1) = RHS_uh + ( A_uh - D_uh)%*%x_i + D_uh%*%x_(i+1)
        !     vh_(i+1) = RHS_vh + ( A_vh - D_vh)%*%x_i + D_vh%*%x_(i+1)
        ! until the differences between x_(i+1) and x_(i) are negligable. The RHS terms:
        !     RHS_uh = [shallow_water_solution_uh - A_uh%*%x_last_timestep ]
        !     RHS_vh = [shallow_water_solution_vh - A_vh%*%x_last_timestep ]
        ! do not change in time.
        !
        !
        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
            dlon, dlat, msl_linear,  &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(out) :: offdiagonalAx_uh(:,:), offdiagonalAx_vh(:,:), diagonalA_uh(:,:), diagonalA_vh(:,:)
        logical, intent(in) :: update_UH, update_VH

        integer(ip) :: i, j
        real(dp) :: r_coslat_dlat, r_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
        real(dp) :: r_dlat, r_coslat_dlon_jp1, r_coslat_dlon_j, r_coslat_dlat_jp1, r_coslat_dlat_j, coslat_jp1h

        ! BEWARE: JAGURS USES A DIFFERENT COORDINATE SYSTEM TO SWALS
        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 2, size(uh, 2) - 1

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and cartesian
                ! coordinates.
                r_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
                r_coslat_dlat   = r_coslat_dlon * dlat / dlon
                coslat_jph = distance_bottom_edge(j+1) * dlat / (dlon * distance_left_edge(1))
                coslat_jmh = distance_bottom_edge(j+0) * dlat / (dlon * distance_left_edge(1))

                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)

                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * r_coslat_dlon)

                    ! The diagonal component of A_uh (related to uh(i,j) )
                    diagonalA_uh(i,j) = dispersive_premult * 1.0_dp/r_coslat_dlon * (-2.0_dp)

                    offdiagonalAx_uh(i,j) = &
                        ! First subtract the diagonal_uh%*%x component
                        -diagonalA_uh(i,j) * uh(i,j) + &
                        ! Then include all terms in A_uh%*%x here (easier bookkeeping to write it this way)
                        dispersive_premult * ( &
                            ! 1/(R_coslat) * d/dlon(uh)
                            1.0_dp / r_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
                            ! 1/(R_coslat) * d/dlat( coslat * vh )
                            1.0_dp / r_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) )
                end do
            end do
            !$OMP END DO
        end if

        if(update_VH) then
            !$OMP DO
            do j = 2, size(vh, 2) - 1

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_dlat = distance_left_edge(1)
                r_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
                r_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
                r_coslat_dlat_jp1 = r_coslat_dlon_jp1 * dlat / dlon
                r_coslat_dlat_j   = r_coslat_dlon_j   * dlat / dlon

                coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( distance_left_edge(1)  * dlon )
                coslat_jph  = distance_bottom_edge(j+1) * dlat / ( distance_left_edge(1)  * dlon )
                coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( distance_left_edge(1)  * dlon )

                do i = 2, size(vh, 1) - 1

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)

                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat)

                    ! The diagonal component of A_vh (related to vh(i,j))
                    diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp / r_coslat_dlat_jp1 * coslat_jph + &
                                                               1.0_dp / r_coslat_dlat_j   * coslat_jph )

                    offdiagonalAx_vh(i,j) = &
                        ! First subtract the diagonal_vh%*%x component
                        -diagonalA_vh(i,j) * vh(i,j) + &
                        ! Then include all terms in A_vh%*%x here (easier bookkeeping to write it this way)
                        dispersive_premult * ( &
                            ! uh term
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) + &
                            ! vh term
                            ( 1.0_dp / r_coslat_dlat_jp1 * ( coslat_jp1h * vh(i, j+1) - coslat_jph * vh(i, j  )) - &
                              1.0_dp / r_coslat_dlat_j   * ( coslat_jph  * vh(i, j  ) - coslat_jmh * vh(i, j-1)) ) )
                end do
            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine

    subroutine dispersive_solve_staggered_grid_JACOBI(ds, U, &
            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear)
        ! Use Jacobi iteration to solve the linear dispersive discretization presented in
        ! Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
        ! Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(inout) :: U(:,:,:)
            ! domain%U after an explicit shallow water update
        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver

        integer :: i, j, iter
        real(dp) :: max_err, last_U

        ! Get the explicit part of the dispersive terms
        call linear_dispersive_matmult_JACOBI(&
            ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
            msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
            ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
            ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
            update_UH=.true., update_VH=.true.)

        ! Setup the right-hand-side term, combining the shallow-water solution with the explicit part
        ! of the dispersive term
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds, U)
        !$OMP DO
        do j = 2, size(U,2) - 1
            ds%RHS(:,j,UH) = U(:,j,UH) - (ds%offdiagonalAx(:,j,UH) + ds%diagonalA(:,j,UH) * ds%last_U(:,j,UH))
        end do
        !$OMP END DO NOWAIT
        !$OMP DO
        do j = 2, size(U,2) - 1
            ds%RHS(:,j,VH) = U(:,j,VH) - (ds%offdiagonalAx(:,j,VH) + ds%diagonalA(:,j,VH) * ds%last_U(:,j,VH))
        end do
        !$OMP END DO

        !$OMP END PARALLEL

        !call ds%store_last_U(U)

        jacobi_iter: do iter = 1, ds%max_iter

            ! Jacobi iteration
            ds%last_iter = iter

            !
            !UH update
            !
            call linear_dispersive_matmult_JACOBI(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.true., update_VH=.false.)

            max_err = 0.0_dp
            !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            do j = 2, size(U,2) - 1
                do i = 2, size(U, 1) - 1
                    last_U = U(i,j,UH)
                    U(i,j,UH) = (ds%RHS(i,j,UH) + ds%offdiagonalAx(i,j,UH))/(1.0_dp - ds%diagonalA(i,j,UH))
                    U(i,j,UH) = U(i,j,UH) + jacobi_overrelax*(U(i,j,UH)-last_U)

                    ! Record the max abs_uh_difference/depth, reducing to abs_uh_difference in depths < 1m
                    max_err = max( max_err, &
                        (abs(U(i,j,UH) - last_U)/&
                         max(msl_linear - 0.5_dp * (U(i+1,j,ELV) + U(i,j,ELV)), 1.0_dp)) )
                end do
            end do
            !$OMP END PARALLEL DO

            !print*, 'Uchange : ', maxval(abs(U(:,:,UH) - ds%last_U(:,:,UH)))

            !
            ! VH update
            !
            call linear_dispersive_matmult_JACOBI(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.false., update_VH=.true.)

            !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            do j = 2, size(U,2) - 1
                do i = 2, size(U, 1) - 1
                    last_U = U(i,j,VH)
                    U(i,j,VH) = (ds%RHS(i,j,VH) + ds%offdiagonalAx(i,j,VH))/(1.0_dp - ds%diagonalA(i,j,VH))
                    U(i,j,VH) = U(i,j,VH) + jacobi_overrelax*(U(i,j,VH)-last_U)
                    ! Record the max abs_vh_difference/depth, reducing to abs_vh_difference in depths < 1m
                    max_err = max( max_err, &
                        (abs(U(i, j,VH) - last_U)/&
                         max(msl_linear - 0.5_dp * (U(i,j,ELV) + U(i,j+1,ELV)), 1.0_dp)) )
                end do
            end do
            !$OMP END PARALLEL DO

            !print*, '      err', max_err, '; iter ', iter

            ! Check for tolerance
            if(max_err < ds%tol) exit jacobi_iter

        end do jacobi_iter

        !print*, 'Jacobi iter: ', ds%last_iter, max_err!, maxval(abs(U(:,:,UH:VH) - ds%last_U(:,:,UH:VH)))
        if(ds%last_iter == ds%max_iter) then
            write(log_output_unit, *) 'Jacobi iteration hit max iterations (', ds%max_iter, ') with error ', max_err
        end if

    end subroutine


#ifdef DISPERSIVE_TRIDIAG
    subroutine linear_dispersive_sweep_TRIDIAG(elev, uh, vh, &
            RHS_uh, RHS_vh, &
            msl_linear, distance_bottom_edge, distance_left_edge, &
            dlon, dlat, update_UH, update_VH)

        real(dp), intent(in) :: elev(:,:), RHS_uh(:,:), RHS_vh(:,:), &
            dlon, dlat, msl_linear,  &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(inout) :: uh(:,:), vh(:,:)
        logical, intent(in) :: update_UH, update_VH

        integer(ip) :: i, j
        integer :: N, INFO
        real(dp) :: r_coslat_dlat, r_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
        real(dp) :: r_dlat, r_coslat_dlon_jp1, r_coslat_dlon_j, r_coslat_dlat_jp1, r_coslat_dlat_j, coslat_jp1h
        ! Arrays for tridagonal solves
        real(dp) :: diagonalA_uh(size(uh, 1)), lowerA_uh(size(uh, 1)), upperA_uh(size(uh, 1)), b_uh(size(uh,1))
        real(dp) :: diagonalA_vh(size(vh, 2)), lowerA_vh(size(vh, 2)), upperA_vh(size(vh, 2)), b_vh(size(vh,2))

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, RHS_uh, RHS_vh, msl_linear, &
        !$OMP                                  dlat, dlon, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 2, size(uh, 2) - 1

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and cartesian
                ! coordinates.
                r_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
                r_coslat_dlat   = r_coslat_dlon * dlat / dlon
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
                b_uh(1) = uh(1,j)
                b_uh(N) = uh(N,j)

                do i = 2, size(uh, 1) - 1
                    ! UH dispersive term

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)

                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * r_coslat_dlon)

                    ! Terms related to the uh derivative (with 1.0 --> time update)
                    diagonalA_uh(i) = 1.0_dp + 2.0_dp * dispersive_premult * 1.0_dp/r_coslat_dlon
                    lowerA_uh(i) = -dispersive_premult * 1.0_dp/r_coslat_dlon
                    upperA_uh(i) = -dispersive_premult * 1.0_dp/r_coslat_dlon

                    b_uh(i) = RHS_uh(i,j) + &
                        dispersive_premult * ( &
                            ! 1/(R_coslat) * d/dlat( coslat * vh )
                            1.0_dp / r_coslat_dlat * ( (coslat_jph * vh(i+1,j) - coslat_jmh * vh(i+1,j-1)) - &
                                                       (coslat_jph * vh(i  ,j) - coslat_jmh * vh(i  ,j-1)) ) )
                end do
                ! Compute the solution with lapack's tridiagonal solver
#ifdef REALFLOAT
                call sgtsv(N, 1, lowerA_uh(2:N), diagonalA_uh, upperA_uh(1:N-1), b_uh, N, INFO)
#else
                call dgtsv(N, 1, lowerA_uh(2:N), diagonalA_uh, upperA_uh(1:N-1), b_uh, N, INFO)
#endif
                !uh(:,j) = b_uh
                uh(:,j) = b_uh + (tridiagonal_overrelax) * (b_uh - uh(:,j))

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
                b_vh(1) = vh(i,1)
                b_vh(N) = vh(i,N)

                do j = 2, size(vh, 2) - 1

                    ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                    ! cartesian coordinates.
                    r_dlat = distance_left_edge(1)
                    r_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
                    r_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
                    r_coslat_dlat_jp1 = r_coslat_dlon_jp1 * dlat / dlon
                    r_coslat_dlat_j   = r_coslat_dlon_j   * dlat / dlon

                    coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( distance_left_edge(1)  * dlon )
                    coslat_jph  = distance_bottom_edge(j+1) * dlat / ( distance_left_edge(1)  * dlon )
                    coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( distance_left_edge(1)  * dlon )

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)

                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat)

                    ! Terms related to the vh derivative (with 1.0 --> time update)
                    diagonalA_vh(j) = 1.0_dp + dispersive_premult * (1.0_dp / r_coslat_dlat_jp1 * coslat_jph + &
                                                                     1.0_dp / r_coslat_dlat_j   * coslat_jph )
                    lowerA_vh(j) = -dispersive_premult * 1.0_dp / r_coslat_dlat_j   * coslat_jmh
                    upperA_vh(j) = -dispersive_premult * 1.0_dp / r_coslat_dlat_jp1 * coslat_jp1h

                    b_vh(j) = RHS_vh(i,j) + &
                        dispersive_premult * ( &
                            ! uh term
                            ( 1.0_dp / R_coslat_dlon_jp1 * (uh(i  , j+1) - uh(i-1, j+1)) - &
                              1.0_dp / R_coslat_dlon_j   * (uh(i  , j  ) - uh(i-1, j  )) ) )
                end do
                ! Compute the solution with lapack's tridiagonal solver
#ifdef REALFLOAT
                call sgtsv(N, 1, lowerA_vh(2:N), diagonalA_vh, upperA_vh(1:N-1), b_vh, N, INFO)
#else
                call dgtsv(N, 1, lowerA_vh(2:N), diagonalA_vh, upperA_vh(1:N-1), b_vh, N, INFO)
#endif
                !vh(i,:) = b_vh
                vh(i,:) = b_vh + (tridiagonal_overrelax) * (b_vh - vh(i,:))

            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine

    subroutine dispersive_solve_staggered_grid_TRIDIAG(ds, U, &
            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear)
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(inout) :: U(:,:,:)
            ! domain%U after an explicit shallow water update
        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver

        integer :: i, j, iter
        real(dp) :: max_err, last_U

        ! Get the explicit part of the dispersive terms
        call linear_dispersive_matmult_JACOBI(&
            ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
            msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
            ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
            ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
            update_UH=.true., update_VH=.true.)

        ! Setup the right-hand-side term, combining the shallow-water solution with the explicit part
        ! of the dispersive term
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds, U)
        !$OMP DO
        do j = 2, size(U,2) - 1
            ds%RHS(:,j,UH) = U(:,j,UH) - (ds%offdiagonalAx(:,j,UH) + ds%diagonalA(:,j,UH) * ds%last_U(:,j,UH))
        end do
        !$OMP END DO NOWAIT
        !$OMP DO
        do j = 2, size(U,2) - 1
            ds%RHS(:,j,VH) = U(:,j,VH) - (ds%offdiagonalAx(:,j,VH) + ds%diagonalA(:,j,VH) * ds%last_U(:,j,VH))
        end do
        !$OMP END DO

        !$OMP END PARALLEL

        iter_loop: do iter = 1, ds%max_iter

            ds%last_iter = iter

            call ds%store_last_U(U)

            !
            !UH and VH update
            !
            call linear_dispersive_sweep_TRIDIAG(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                ds%RHS(:,:,UH), ds%RHS(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                update_UH=.true., update_VH=.true.)

            max_err = 0.0_dp
            !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            !$OMP DO
            do j = 2, size(U,2) - 1
                do i = 2, size(U, 1) - 1
                    ! Record the max abs_uh_difference/depth, reducing to abs_uh_difference in depths < 1m
                    max_err = max( max_err, &
                        (abs(U(i,j,UH) - ds%last_U(i,j,UH) )/&
                         max(msl_linear - 0.5_dp * (U(i+1,j,ELV) + U(i,j,ELV)), 1.0_dp)) )
                end do
            end do
            !$OMP END DO NOWAIT

            !$OMP DO
            do j = 2, size(U,2) - 1
                do i = 2, size(U, 1) - 1
                    ! Record the max abs_vh_difference/depth, reducing to abs_vh_difference in depths < 1m
                    max_err = max( max_err, &
                        (abs(U(i, j,VH) - ds%last_U(i,j,VH) )/&
                         max(msl_linear - 0.5_dp * (U(i,j,ELV) + U(i,j+1,ELV)), 1.0_dp)) )
                end do
            end do
            !$OMP END DO
            !$OMP END PARALLEL

            !print*, '       ', max_err

            ! Check for tolerance
            if(max_err < ds%tol) exit iter_loop

        end do iter_loop

        !print*, 'Tridiagonal iter: ', ds%last_iter, max_err
        if(ds%last_iter == ds%max_iter) then
            write(log_output_unit, *) 'Tridiagonal iteration hit max iterations (', ds%max_iter, ') with error ', max_err
        end if

    end subroutine
#endif
end module
