! Multigrid ideas might help
! https://petsc.org/release/src/ksp/ksp/tutorials/ex22f.F90.html

! Use definition below to use jacobi iteration
!#define DISPERSIVE_JACOBI

! Use definition below to use spatially variable jacobi iteration. This skips
! iterations in areas with errors < threshold, with some overhead for the adaptivity.
#define DISPERSIVE_JACOBI_ADAPTIVE

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
    !real(dp), parameter :: td1 = 0.1_dp, td2 = 0.09999_dp ! Linear taper of dispersive terms at depth
    ! real(dp), parameter :: td1 = 0._dp, td2 = -1.0_dp ! TURN OFF LINEAR TAPER of dispersive terms at depth

    type dispersive_solver_type

        logical :: is_staggered_grid 
            !! Different solvers for staggered vs cell-centred grids
        real(dp), allocatable :: RHS(:,:,:), last_U(:,:,:), offdiagonalAx(:,:,:), diagonalA(:,:,:)
            !! Work arrays

#ifdef DISPERSIVE_JACOBI_ADAPTIVE        
        ! These help with spatially variable jacobi iteration
        real(dp), allocatable :: local_err(:,:)
        logical, allocatable :: needs_update(:,:)
        integer(ip), allocatable :: i0(:), i1(:)
        real(dp) :: local_tol_factor = 0.5_dp
                !! To stop iterating locally we require the change in the solution to be less
                !! than "local_tol_factor * (the solver tolerance)", i.e., substantially less than the 
                !! tolerance. 
#endif

#ifdef REALFLOAT
        real(dp) :: tol = 1.0e-06_dp !! Solver tolerance -- iteration until this accuracy is reached
#else
        real(dp) :: tol = 1.0e-08_dp !! Solver tolerance -- iteration until this accuracy is reached
#endif
        integer(ip) :: max_iter = 2000 !! Maximum number of iterations in a single call to ds%solve
        integer(ip) :: last_iter = 0 !! Iteration count at convergence
        real(dp) :: max_err = -HUGE(1.0_dp) !! Store the max error on the last iteration

        real(dp) :: td1 = 0.0_dp, td2 = -1.0_dp 
            !! Depths (below MSL) used to linearly taper-off dispersion. 
            !! Default corresponds to no tapering. 
            !! Set td1 > td2 to linearly reduce dispersive terms to zero between depths of td1 --> td2

        contains
        procedure :: store_last_U => store_last_U
        procedure :: setup => setup_dispersive_solver
        procedure :: solve => linear_dispersive_solve

#ifdef DISPERSIVE_JACOBI
        procedure :: solve_staggered_grid => linear_dispersive_solve_staggered_grid_JACOBI
#endif

#ifdef DISPERSIVE_JACOBI_ADAPTIVE
        procedure :: solve_staggered_grid => linear_dispersive_solve_staggered_grid_JACOBI_ADAPTIVE
        procedure :: solve_cellcentred_grid => linear_dispersive_solve_cellcentred_JACOBI_ADAPTIVE
#endif

#ifdef DISPERSIVE_TRIDIAG
        procedure :: solve_staggered_grid => linear_dispersive_solve_staggered_grid_TRIDIAG
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

    subroutine setup_dispersive_solver(ds, nx, ny, is_staggered_grid)
        !! Allocate workspace
        class(dispersive_solver_type), intent(inout) :: ds
        integer(ip), intent(in) :: nx, ny
        logical, intent(in) :: is_staggered_grid

        integer(ip) :: i, j

        ds%is_staggered_grid = is_staggered_grid

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
            ds%diagonalA(:,j,UH:VH) = 1.0_dp
        end do
        !$OMP END PARALLEL DO

#ifdef DISPERSIVE_JACOBI_ADAPTIVE
        if(allocated(ds%local_err)) deallocate(ds%local_err)
        if(allocated(ds%needs_update)) deallocate(ds%needs_update)
        if(allocated(ds%i0)) deallocate(ds%i0)
        if(allocated(ds%i1)) deallocate(ds%i1)

        allocate(ds%local_err(nx, ny), ds%needs_update(nx, ny), ds%i0(ny), ds%i1(ny))
        !$OMP PARALLEL DO
        do j = 1, ny
            ds%local_err(:,j) = 0.0_dp
            ds%needs_update(:,j) = .true.
            ds%i0(j) = 1
            ds%i1(j) = nx
        end do
        !$OMP END PARALLEL DO
#endif

    end subroutine

    subroutine linear_dispersive_solve(ds, U, &
            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear, rhs_is_up_to_date)
        !! Convenience wrapper which gives a consistent interface for both staggered and cell-centred solvers
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(inout) :: U(:,:,:)
            ! domain%U after an explicit shallow water update
        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver
        logical, optional, intent(in) :: rhs_is_up_to_date

        if(ds%is_staggered_grid) then
            ! Staggered grid dispersive solve
            call ds%solve_staggered_grid(U, dlon, dlat, &
                distance_bottom_edge, distance_left_edge, &
                msl_linear, rhs_is_up_to_date)
        else
            ! Cell-centred dispersive solve
            call ds%solve_cellcentred_grid(U, dlon, dlat, &
                distance_bottom_edge, distance_left_edge, &
                msl_linear, rhs_is_up_to_date)
        end if

    end subroutine

    !
    !
    !
    subroutine linear_dispersive_staggered_matmult(elev, uh, vh, msl_linear, distance_bottom_edge, distance_left_edge, &
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

    !subroutine linear_dispersive_cellcentred_matmult(elev, uh, vh, msl_linear, distance_bottom_edge, distance_left_edge, &
    !        dlon, dlat, solution_uh, solution_vh, update_UH, update_VH)
    !    !
    !    ! Towards solving the linear dispersive equation as per:
    !    !   Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
    !    !       Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
    !    ! Or
    !    !   Baba, T.; Allgeyer, S.; Hossen, J.; Cummins, P. R.; Tsushima, H.; Imai, K.; Yamashita, K. & Kato, T. Accurate numerical
    !    !   simulation of the far-field tsunami caused by the 2011 Tohoku earthquake, including the effects of Boussinesq dispersion,
    !    !   seawater density stratification, elastic loading, and gravitational potential change Ocean Modelling, Elsevier BV, 2017,
    !    !   111, 46–54
    !    !
    !    ! BACKGROUND:
    !    !   This routine is a cell-centred varient of subroutine "linear_dispersive_staggered_matmult". See documentation there
    !    !   for more information
    !    !
    !    real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
    !        dlon, dlat, msl_linear,  &
    !        distance_bottom_edge(:), distance_left_edge(:)
    !    real(dp), intent(out) :: solution_uh(:,:), solution_vh(:,:)
    !    logical, intent(in) :: update_UH, update_VH

    !    integer(ip) :: i, j
    !    real(dp) :: r_coslat_dlat, r_coslat_dlon, coslat_jp1, coslat_j, coslat_jm1, d_i_j, dispersive_premult
    !    real(dp) :: r_dlat, r_coslat_dlon_jph, r_coslat_dlon_jmh, r_coslat_dlat_jph, r_coslat_dlat_jmh

    !    !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, solution_uh, solution_vh, &
    !    !$OMP                                  dlat, dlon, &
    !    !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

    !    if(update_UH) then
    !        !$OMP DO
    !        do j = 1, size(uh, 2)

    !            if(j == 1 .or. j == size(uh, 2)) then
    !                solution_uh(:,j) = uh(:,j)
    !                cycle
    !            end if

    !            ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
    !            ! cartesian coordinates.
    !            r_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  ))
    !            r_coslat_dlat   = r_coslat_dlon * dlat / dlon

    !            ! NOTE regarding array bounds. 
    !            ! - We know j > 1 and j < size(uh,2). 
    !            ! - Also distance_bottom_edge has size == size(uh,2)+1
    !            coslat_jp1 = 0.5_dp * (distance_bottom_edge(j+2) + distance_bottom_edge(j+1)) * &
    !                dlat / (dlon * distance_left_edge(1))
    !            coslat_j   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j+0)) * &
    !                dlat / (dlon * distance_left_edge(1))
    !            coslat_jm1 = 0.5_dp * (distance_bottom_edge(j+0) + distance_bottom_edge(j-1)) * &
    !                dlat / (dlon * distance_left_edge(1))

    !            solution_uh(1,j) = uh(1,j)
    !            solution_uh(size(uh,1),j) = uh(size(uh,1),j)
    !            do i = 2, size(uh, 1) - 1
    !                ! UH dispersive term

    !                ! Depth at location of uh(i,j)
    !                !d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
    !                !    elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
    !                d_i_j = merge(msl_linear - elev(i, j), 0.0_dp, elev(i,j) < msl_linear)

    !                !dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * r_coslat_dlon)
    !                dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_coslat_dlon)

    !                solution_uh(i,j) = uh(i,j) - &
    !                    ! Include all terms in A_uh%*%x here
    !                    ! h0^2 / (3 R_coslat) * d/dlon
    !                    dispersive_premult * ( &
    !                        ! lon-diff { 1/(R_coslat) * d/dlon(uh) }
    !                        1.0_dp / r_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
    !                        ! lon-diff { 1/(R_coslat) * d/dlat( coslat * vh ) }
    !                        1.0_dp / r_coslat_dlat * ( &
    !                            (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
    !                            (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp)  &
    !                            )
    !            end do
    !        end do
    !        !$OMP END DO NOWAIT
    !    end if

    !    if(update_VH) then
    !        !$OMP DO
    !        do j = 2, size(vh, 2) - 1

    !            if(j == 1 .or. j == size(vh, 2)) then
    !                solution_vh(:,j) = vh(:,j)
    !                cycle
    !            end if

    !            ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
    !            ! cartesian coordinates.
    !            r_dlat = distance_left_edge(1)
    !            r_coslat_dlon_jph = distance_bottom_edge(j+1)
    !            r_coslat_dlon_jmh = distance_bottom_edge(j+0)
    !            r_coslat_dlat_jph = r_coslat_dlon_jph * dlat / dlon
    !            r_coslat_dlat_jmh = r_coslat_dlon_jmh * dlat / dlon

    !            coslat_jp1 = 0.5_dp*(distance_bottom_edge(j+2)+distance_bottom_edge(j+1)) * dlat / (distance_left_edge(1) * dlon)
    !            coslat_j   = 0.5_dp*(distance_bottom_edge(j+1)+distance_bottom_edge(j+0)) * dlat / (distance_left_edge(1) * dlon)
    !            coslat_jm1 = 0.5_dp*(distance_bottom_edge(j+0)+distance_bottom_edge(j-1)) * dlat / (distance_left_edge(1) * dlon)

    !            solution_vh(1,j) = vh(1,j)
    !            solution_vh(size(vh,1), j) = vh(size(vh,1), j)
    !            do i = 2, size(vh, 1) - 1

    !                ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
    !                !d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
    !                !    elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
    !                d_i_j = merge(msl_linear - elev(i,j), 0.0_dp, elev(i,j) < msl_linear)

    !                dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_dlat)

    !                solution_vh(i,j) = vh(i,j) - &
    !                    ! Include all terms in A_vh%*%x here
    !                    ! h0**2 / (3 * R) * d/dlat * (
    !                    dispersive_premult * ( &
    !                        ! uh term
    !                        ! 1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of centred differences at j+1 and j.
    !                        ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
    !                        ! 1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of centred differences at j and j-1
    !                          1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) + &
    !                        ! vh term
    !                        ! 1/(R coslat ) d(coslat vh)/dlat @ i, j+1/2
    !                        ( 1.0_dp / r_coslat_dlat_jph * ( coslat_jp1 * vh(i, j+1) - coslat_j   * vh(i, j  )) - &
    !                        ! 1/(R coslat ) d(coslat vh)/dlat @ i, j-1/2
    !                          1.0_dp / r_coslat_dlat_jmh * ( coslat_j   * vh(i, j  ) - coslat_jm1 * vh(i, j-1)) ) )
    !            end do
    !        end do
    !        !$OMP END DO
    !    end if

    !    !$OMP END PARALLEL

    !end subroutine

    !
    ! JACOBI ITERATION
    !

    subroutine linear_dispersive_staggered_matmult_JACOBI(elev, uh, vh, msl_linear, &
            distance_bottom_edge, distance_left_edge, &
            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, &
            update_UH, update_VH)
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

        ! BEWARE: THE JAGURS PAPER USES A DIFFERENT COORDINATE SYSTEM TO SWALS
        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH)

        if(update_UH) then
            !$OMP DO
            do j = 2, size(uh, 2) - 1

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
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

                coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( r_dlat  * dlon )
                coslat_jph  = distance_bottom_edge(j+1) * dlat / ( r_dlat  * dlon )
                coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( r_dlat  * dlon )

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

    subroutine linear_dispersive_solve_staggered_grid_JACOBI(ds, U, &
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
        call linear_dispersive_staggered_matmult_JACOBI(&
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
            call linear_dispersive_staggered_matmult_JACOBI(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
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
            call linear_dispersive_staggered_matmult_JACOBI(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
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

        !!print*, 'Jacobi iter: ', ds%last_iter, max_err!, maxval(abs(U(:,:,UH:VH) - ds%last_U(:,:,UH:VH)))
        !if(ds%last_iter == ds%max_iter) then
        !    write(log_output_unit, *) 'Jacobi iteration hit max iterations (', ds%max_iter, ') with error ', max_err
        !end if

    end subroutine


#ifdef DISPERSIVE_JACOBI_ADAPTIVE

    subroutine linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(elev, uh, vh, msl_linear, &
            distance_bottom_edge, distance_left_edge, &
            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, &
            update_UH, update_VH, needs_update, i0, i1, td1, td2)
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
        logical, intent(in) :: update_UH, update_VH, needs_update(:,:)
        integer(ip), intent(in) :: i0(:), i1(:)
        real(dp), intent(in) :: td1, td2 ! "Depths below msl" used to taper off dispersion

        integer(ip) :: i, j
        real(dp) :: r_coslat_dlat, r_coslat_dlon, coslat_jph, coslat_jmh, d_iph_j, d_i_jph, dispersive_premult
        real(dp) :: r_dlat, r_coslat_dlon_jp1, r_coslat_dlon_j, r_coslat_dlat_jp1, r_coslat_dlat_j, coslat_jp1h

        ! BEWARE: THE JAGURS PAPER USES A DIFFERENT COORDINATE SYSTEM TO SWALS
        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, needs_update, i0, i1, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH, td1, td2)

        if(update_UH) then
            !$OMP DO SCHEDULE(DYNAMIC, 10)
            !!$OMP DO SCHEDULE(STATIC)
            do j = 2, size(uh, 2) - 1
                !if(.not. any(needs_update(:,j))) cycle

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  )) ! x-cell-distance at uh(i,j)
                r_coslat_dlat   = r_coslat_dlon * dlat / dlon
                coslat_jph = distance_bottom_edge(j+1) * dlat / (dlon * distance_left_edge(1))
                coslat_jmh = distance_bottom_edge(j+0) * dlat / (dlon * distance_left_edge(1))

                !do i = 2, size(uh, 1) - 1
                do i = i0(j), i1(j)
                    ! UH dispersive term
                    if(.not. needs_update(i,j)) cycle

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)

                    dispersive_premult = d_iph_j*d_iph_j / (3.0_dp * r_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_iph_j - td2)/(td1 - td2)))

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
            !$OMP DO SCHEDULE(DYNAMIC, 10)
            !!$OMP DO SCHEDULE(STATIC)
            do j = 2, size(vh, 2) - 1
                !if(.not. any(needs_update(:,j))) cycle

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_dlat = distance_left_edge(1)
                r_coslat_dlon_jp1 = 0.5_dp * ( distance_bottom_edge(j+2) + distance_bottom_edge(j+1) )
                r_coslat_dlon_j   = 0.5_dp * ( distance_bottom_edge(j+1) + distance_bottom_edge(j+0) )
                r_coslat_dlat_jp1 = r_coslat_dlon_jp1 * dlat / dlon
                r_coslat_dlat_j   = r_coslat_dlon_j   * dlat / dlon

                coslat_jp1h = distance_bottom_edge(j+2) * dlat / ( r_dlat  * dlon )
                coslat_jph  = distance_bottom_edge(j+1) * dlat / ( r_dlat  * dlon )
                coslat_jmh  = distance_bottom_edge(j+0) * dlat / ( r_dlat  * dlon )

                !do i = 2, size(vh, 1) - 1
                do i = i0(j), i1(j)
                    if(.not. needs_update(i,j)) cycle

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                        elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)

                    dispersive_premult = d_i_jph*d_i_jph / (3.0_dp * r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_jph - td2)/(td1 - td2)))


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

    subroutine linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(elev, uh, vh, msl_linear, &
            distance_bottom_edge, distance_left_edge, &
            dlon, dlat, offdiagonalAx_uh, offdiagonalAx_vh, diagonalA_uh, diagonalA_vh, &
            update_UH, update_VH, needs_update, i0, i1, td1, td2)
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
        !   This is a cell-centred version of linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE
        !   See that routine for documentation
        !
        !
        real(dp), intent(in) :: elev(:,:), uh(:,:), vh(:,:), &
            dlon, dlat, msl_linear,  &
            distance_bottom_edge(:), distance_left_edge(:)
        real(dp), intent(out) :: offdiagonalAx_uh(:,:), offdiagonalAx_vh(:,:), diagonalA_uh(:,:), diagonalA_vh(:,:)
        logical, intent(in) :: update_UH, update_VH, needs_update(:,:)
        integer(ip), intent(in) :: i0(:), i1(:)
        real(dp), intent(in) :: td1, td2 ! "Depth below MSL" used to linear taper off dispersive terms

        integer(ip) :: i, j
        real(dp) :: r_coslat_dlat, r_coslat_dlon, coslat_jp1, coslat_j, coslat_jm1, d_i_j, dispersive_premult
        real(dp) :: r_dlat, r_coslat_dlon_jph, r_coslat_dlon_jmh, r_coslat_dlat_jph, r_coslat_dlat_jmh

        ! BEWARE: THE JAGURS PAPER USES A DIFFERENT COORDINATE SYSTEM TO SWALS
        ! They use co-latitude -- so sin(lat) becomes cos(lat) in SWALS coordinates.

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(uh, vh, elev, msl_linear, offdiagonalAx_uh, offdiagonalAx_vh, &
        !$OMP                                  diagonalA_uh, diagonalA_vh, dlat, dlon, needs_update, i0, i1, &
        !$OMP                                  distance_left_edge, distance_bottom_edge, update_UH, update_VH, td1, td2)

        if(update_UH) then
            !$OMP DO SCHEDULE(DYNAMIC, 10)
            !!$OMP DO SCHEDULE(STATIC)
            do j = 2, size(uh, 2) - 1
                !if(.not. any(needs_update(:,j))) cycle

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_coslat_dlon   = 0.5_dp * (distance_bottom_edge(j+1) + distance_bottom_edge(j  ))
                r_coslat_dlat   = r_coslat_dlon * dlat / dlon

                ! NOTE regarding array bounds. 
                ! - We know j > 1 and j < size(uh,2). 
                ! - Also distance_bottom_edge has size == size(uh,2)+1
                coslat_jp1 = 0.5_dp * (distance_bottom_edge(j+2) + distance_bottom_edge(j+1)) * &
                    dlat / (dlon * distance_left_edge(1))
                coslat_jm1 = 0.5_dp * (distance_bottom_edge(j+0) + distance_bottom_edge(j-1)) * &
                    dlat / (dlon * distance_left_edge(1))

                !do i = 2, size(uh, 1) - 1
                do i = i0(j), i1(j)
                    ! UH dispersive term
                    if(.not. needs_update(i,j)) cycle

                    ! Depth at location of uh(i,j) -- or zero if either neighbour is dry
                    !d_iph_j = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i+1,j)), 0.0_dp, &
                    !    elev(i,j) < msl_linear .and. elev(i+1, j) < msl_linear)
                    d_i_j = merge(msl_linear - elev(i, j), 0.0_dp, elev(i,j) < msl_linear)

                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_coslat_dlon) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! The diagonal component of A_uh (related to uh(i,j) )
                    diagonalA_uh(i,j) = dispersive_premult * 1.0_dp/r_coslat_dlon * (-2.0_dp)

                    offdiagonalAx_uh(i,j) = &
                        ! First subtract the diagonal_uh%*%x component
                        -diagonalA_uh(i,j) * uh(i,j) + &
                        ! Then include all terms in A_uh%*%x here (easier bookkeeping to write it this way)
                        ! h0^2 / (3 R_coslat) * d/dlon
                        dispersive_premult * ( &
                            ! lon-diff of { 1/(R_coslat) * d/dlon(uh) }
                            1.0_dp / r_coslat_dlon * ( uh(i+1,j) - 2.0_dp * uh(i,j) + uh(i-1,j) ) + &
                            ! lon-diff of { 1/(R_coslat) * d/dlat( coslat * vh ) }
                            1.0_dp / r_coslat_dlat * ( &
                                (coslat_jp1 * vh(i+1,j+1) - coslat_jm1 * vh(i+1,j-1))*0.25_dp - &
                                (coslat_jp1 * vh(i-1,j+1) - coslat_jm1 * vh(i-1,j-1))*0.25_dp)  &
                                )

                end do
            end do
            !$OMP END DO
        end if

        if(update_VH) then
            !$OMP DO SCHEDULE(DYNAMIC, 10)
            !!$OMP DO SCHEDULE(STATIC)
            do j = 2, size(vh, 2) - 1
                !if(.not. any(needs_update(:,j))) cycle

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
                r_dlat = distance_left_edge(1)
                r_coslat_dlon_jph = distance_bottom_edge(j+1)
                r_coslat_dlon_jmh = distance_bottom_edge(j+0)
                r_coslat_dlat_jph = r_coslat_dlon_jph * dlat / dlon
                r_coslat_dlat_jmh = r_coslat_dlon_jmh * dlat / dlon

                coslat_jp1 = 0.5_dp*(distance_bottom_edge(j+2)+distance_bottom_edge(j+1)) * dlat / (distance_left_edge(1) * dlon)
                coslat_j   = 0.5_dp*(distance_bottom_edge(j+1)+distance_bottom_edge(j+0)) * dlat / (distance_left_edge(1) * dlon)
                coslat_jm1 = 0.5_dp*(distance_bottom_edge(j+0)+distance_bottom_edge(j-1)) * dlat / (distance_left_edge(1) * dlon)


                !do i = 2, size(vh, 1) - 1
                do i = i0(j), i1(j)
                    if(.not. needs_update(i,j)) cycle

                    ! Depth at location of vh(i,j) -- or zero if either neighbour is dry
                    !d_i_jph = merge(0.5_dp * (msl_linear - elev(i,j) + msl_linear - elev(i,j+1)), 0.0_dp, &
                    !    elev(i,j) < msl_linear .and. elev(i, j+1) < msl_linear)
                    d_i_j = merge(msl_linear - elev(i,j), 0.0_dp, elev(i,j) < msl_linear)

                    dispersive_premult = d_i_j*d_i_j / (3.0_dp * r_dlat) * &
                        ! Linear taper
                        min(1.0_dp, max(0.0_dp, (d_i_j - td2)/(td1 - td2)))

                    ! The diagonal component of A_vh (related to vh(i,j))
                    !diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp / r_coslat_dlat_jp1 * coslat_jph + &
                    !                                           1.0_dp / r_coslat_dlat_j   * coslat_jph )
                    diagonalA_vh(i,j) = -dispersive_premult * (1.0_dp/r_coslat_dlat_jph * coslat_j  + &
                                                               1.0_dp/r_coslat_dlat_jmh * coslat_j )

                    offdiagonalAx_vh(i,j) = &
                        ! First subtract the diagonal_vh%*%x component
                        -diagonalA_vh(i,j) * vh(i,j) + &
                        ! Then include all terms in A_vh%*%x here (easier bookkeeping to write it this way)
                        ! h0**2 / (3 * R) * d/dlat * (
                        dispersive_premult * ( &
                            ! uh term
                            ! 1/(R coslat) d(uh)/dlon @ i, j+1/2, by mean of central differences at j+1 and j.
                            ( 1.0_dp / R_coslat_dlon_jph * (uh(i+1, j+1) - uh(i-1,j+1) + uh(i+1,j) - uh(i-1,j))*0.25_dp - &
                            ! 1/(R coslat) d(uh)/dlon @ i, j-1/2, by mean of central differences at j and j-1
                              1.0_dp / R_coslat_dlon_jmh * (uh(i+1, j-1) - uh(i-1,j-1) + uh(i+1,j) - uh(i-1,j))*0.25_dp ) + &
                            ! vh term
                            ! 1/(R coslat ) d(coslat vh)/dlat @ i, j+1/2
                            ( 1.0_dp / r_coslat_dlat_jph * ( coslat_jp1 * vh(i, j+1) - coslat_j   * vh(i, j  )) - &
                            ! 1/(R coslat ) d(coslat vh)/dlat @ i, j-1/2
                              1.0_dp / r_coslat_dlat_jmh * ( coslat_j   * vh(i, j  ) - coslat_jm1 * vh(i, j-1)) ) )
                end do
            end do
            !$OMP END DO
        end if

        !$OMP END PARALLEL

    end subroutine

    subroutine linear_dispersive_solve_staggered_grid_JACOBI_ADAPTIVE(ds, U, &
            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear, rhs_is_up_to_date)
        ! Use Jacobi iteration to solve the linear dispersive discretization presented in
        ! Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
        ! Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
        !
        ! This version is adaptive -- it skips computations where the error is (locally) small enough.
        !
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(inout) :: U(:,:,:)
            ! domain%U after an explicit shallow water update
        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver
        logical, optional, intent(in) :: rhs_is_up_to_date

        integer :: i, j, iter, i0, j0, i1, j1, ii0, ii1, l0, l1
        real(dp) :: max_err, last_U
        integer, parameter :: min_jacobi_iterations = 1, bw = 2
        logical :: allowed_to_exit, rhs_ready

        rhs_ready = .FALSE.
        if(present(rhs_is_up_to_date)) rhs_ready = rhs_is_up_to_date

        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j = 1, size(U,2)
            ! Ensure that some iteration is performed
            ds%needs_update(:,j) = .true.
            ds%local_err(:,j) = 0.0_dp

            ! i indices to update for each j -- initially look at all cells where we can compute central differences.
            ds%i0(j) = 2
            ds%i1(j) = size(U,1) - 1
        end do
        !$OMP END PARALLEL DO

        if(.not. rhs_ready) then

            ! Get the explicit part of the dispersive terms
            call linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(&
                ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.true., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
                td1 = ds%td1, td2=ds%td2)

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
            !$OMP END DO NOWAIT

            !$OMP END PARALLEL
        end if

        !call ds%store_last_U(U)

        allowed_to_exit = .FALSE. ! FIXME -- do we need this?

        jacobi_iter: do iter = 1, ds%max_iter

            ! Jacobi iteration
            ds%last_iter = iter

            !
            !UH update
            !
            call linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.true., update_VH=.false., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
                td1 = ds%td1, td2 = ds%td2)

            max_err = 0.0_dp
            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            !!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            do j = 2, size(U,2) - 1
                !do i = 2, size(U, 1) - 1
                do i = ds%i0(j), ds%i1(j)
                    if(.not. ds%needs_update(i,j)) cycle

                    last_U = U(i,j,UH)
                    U(i,j,UH) = (ds%RHS(i,j,UH) + ds%offdiagonalAx(i,j,UH))/(1.0_dp - ds%diagonalA(i,j,UH))
                    U(i,j,UH) = U(i,j,UH) + jacobi_overrelax*(U(i,j,UH)-last_U)

                    ! Record the max abs_uh_difference/depth, reducing to abs_uh_difference in depths < 1m
                    ds%local_err(i,j) = ( abs(U(i,j,UH) - last_U)/&
                         max(msl_linear - 0.5_dp * (U(i+1,j,ELV) + U(i,j,ELV)), 1.0_dp) )
                    max_err = max( max_err, ds%local_err(i,j) ) 
                end do
            end do
            !$OMP END PARALLEL DO

            !
            ! VH update
            !
            call linear_dispersive_staggered_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.false., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
                td1 = ds%td1, td2 = ds%td2)

            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            !!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            do j = 2, size(U,2) - 1
                !do i = 2, size(U, 1) - 1
                do i = ds%i0(j), ds%i1(j)
                    if(.not. ds%needs_update(i,j)) cycle

                    last_U = U(i,j,VH)
                    U(i,j,VH) = (ds%RHS(i,j,VH) + ds%offdiagonalAx(i,j,VH))/(1.0_dp - ds%diagonalA(i,j,VH))
                    U(i,j,VH) = U(i,j,VH) + jacobi_overrelax*(U(i,j,VH)-last_U)

                    ! Record the max abs_vh_difference/depth, reducing to abs_vh_difference in depths < 1m
                    ds%local_err(i,j) = max(ds%local_err(i,j), abs(U(i, j,VH) - last_U)/&
                         max(msl_linear - 0.5_dp * (U(i,j,ELV) + U(i,j+1,ELV)), 1.0_dp))
                    max_err = max( max_err, ds%local_err(i,j) )
                end do
            end do
            !$OMP END PARALLEL DO

            ds%max_err = max_err

            ! Check for tolerance
            if(max_err < ds%tol) then
                if(allowed_to_exit) exit jacobi_iter
                ! Once the error is less than the tolerance, we 
                ! start doing full grid iterations again. 
                ! Ideally we still meet the tolerance after one, and ensure 
                ! the tolerance is met everywhere.
                ! 
                allowed_to_exit = .true.
                ds%i0 = 2
                ds%i1 = size(U, 1)-1
                !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds)
                do j = 1, size(U,2)
                    ds%needs_update(:,j) = .true.
                end do
                !$OMP END PARALLEL DO
            end if

            if(iter > (min_jacobi_iterations - 1) .and. (.not. allowed_to_exit)) then
                ! Determine which areas need updates.
                ! We record cells if cells are 'near' a cell that has error > "some fraction of tol" (ds%needs_update(:,:)).
                ! We also record, for each j, a range of i-indices that include all those needing updates.
                !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds)

                !$OMP DO SCHEDULE(DYNAMIC, 10)
                do j = 2, size(U,2) - 1
                    ! Update the 'needs_update' variable

                    ! We don't need to scan every i-index in the j'th row
                    ! Need to look look at range of active cells 'bw' rows above/below the current row,
                    ! as these might have error > tol. Need care to not go out of bounds. 
                    l0 = max(2, j-bw) ! safe j index below
                    l1 = min(size(U,2)-1, j+bw) ! safe j index above

                    ! i indices that we should look at in the j'th row
                    ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
                    ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above

                    do i = ii0 , ii1 
                        ! Check within local window for error > tol
                        i0 = max(1, i-bw)
                        i1 = min(size(U,1), i+bw)
                        j0 = max(1, j-bw)
                        j1 = min(size(U,2), j+bw)
                        ds%needs_update(i, j) = any(ds%local_err(i0:i1,j0:j1) > (ds%tol*ds%local_tol_factor))
                    end do
                end do
                !$OMP END DO

                !$OMP DO SCHEDULE(DYNAMIC, 10)
                do j = 2, size(U, 2) - 1
                    ! Update the 'ds%i0, ds%i1' variables defining the i-range in which
                    ! we search for each j. Could not do this in the previous loop due to
                    ! race condition.

                    ! See the loop above for documentation of the next 4 lines
                    l0 = max(2, j-bw) ! safe j index below
                    l1 = min(size(U,2)-1, j+bw) ! safe j index above
                    ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
                    ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above

                    ! To update in loop -- 
                    ! i0 = min index to update; initialise with a large value
                    ds%i0(j) = size(U, 1)
                    ! i1 = max index to update; initialise with a small value
                    ds%i1(j) = 1
                    do i = ii0, ii1
                        if(ds%needs_update(i,j)) then
                            ! Cell i,j is at least close to a cell that has error above tolerance.
                            ds%i0(j) = min(ds%i0(j), i)
                            ds%i1(j) = max(ds%i1(j), i)
                        end if
                    end do
                end do
                !$OMP END PARALLEL
            end if

        end do jacobi_iter

        !!print*, 'Jacobi iter: ', ds%last_iter, max_err, '; shape(U) = ', shape(U) !count(ds%needs_update)*1.0_dp/(size(U,1)*size(U,2))
        !if(ds%last_iter == ds%max_iter) then
        !    write(log_output_unit, *) 'Jacobi iteration hit max iterations (', ds%max_iter, ') with error ', max_err
        !end if

    end subroutine

    subroutine linear_dispersive_solve_cellcentred_JACOBI_ADAPTIVE(ds, U, &
            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear, rhs_is_up_to_date)
        ! Use Jacobi iteration to solve the linear dispersive discretization presented in
        ! Baba, T.; Takahashi, N.; Kaneda, Y.; Ando, K.; Matsuoka, D. & Kato, T. Parallel Implementation of Dispersive Tsunami Wave
        ! Modeling with a Nesting Algorithm for the 2011 Tohoku Tsunami Pure and Applied Geophysics, 2015, 172, 3455-3472
        !
        ! Cellcentred version of linear_dispersive_solve_staggered_grid_JACOBI_ADAPTIVE
        !
        ! This version is adaptive -- it skips computations where the error is (locally) small enough.
        !
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(inout) :: U(:,:,:)
            ! domain%U after an explicit shallow water update
        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver
        logical, optional, intent(in) :: rhs_is_up_to_date

        integer :: i, j, iter, i0, j0, i1, j1, ii0, ii1, l0, l1
        real(dp) :: max_err, last_U
        integer, parameter :: min_jacobi_iterations = 1, bw = 2
        logical :: allowed_to_exit, rhs_ready

        rhs_ready = .FALSE.
        if(present(rhs_is_up_to_date)) rhs_ready = rhs_is_up_to_date

        !$OMP PARALLEL DO DEFAULT(SHARED)
        do j = 1, size(U,2)
            ! Ensure that some iteration is performed
            ds%needs_update(:,j) = .true.
            ds%local_err(:,j) = 0.0_dp

            ! i indices to update for each j -- initially look at all cells where we can compute central differences.
            ds%i0(j) = 2
            ds%i1(j) = size(U,1) - 1
        end do
        !$OMP END PARALLEL DO

        if(.not. rhs_ready) then

            ! Get the explicit part of the dispersive terms
            call linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(&
                ds%last_U(:,:,ELV), ds%last_U(:,:,UH), ds%last_U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.true., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
                td1 = ds%td1, td2 = ds%td2)

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
            !$OMP END DO NOWAIT

            !$OMP END PARALLEL
        end if
        !call ds%store_last_U(U)

        allowed_to_exit = .FALSE. ! FIXME: do we need this?

        jacobi_iter: do iter = 1, ds%max_iter

            ! Jacobi iteration
            ds%last_iter = iter

            !
            !UH update
            !
            call linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.true., update_VH=.false., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
                td1 = ds%td1, td2 = ds%td2)

            max_err = 0.0_dp
            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            !!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            do j = 2, size(U,2) - 1
                !do i = 2, size(U, 1) - 1
                do i = ds%i0(j), ds%i1(j)
                    if(.not. ds%needs_update(i,j)) cycle

                    last_U = U(i,j,UH)
                    U(i,j,UH) = (ds%RHS(i,j,UH) + ds%offdiagonalAx(i,j,UH))/(1.0_dp - ds%diagonalA(i,j,UH))
                    U(i,j,UH) = U(i,j,UH) + jacobi_overrelax*(U(i,j,UH)-last_U)

                    ! Record the max abs_uh_difference/depth, reducing to abs_uh_difference in depths < 1m
                    ds%local_err(i,j) = ( abs(U(i,j,UH) - last_U)/&
                         max(msl_linear - U(i,j,ELV), 1.0_dp) )
                    max_err = max( max_err, ds%local_err(i,j) ) 
                end do
            end do
            !$OMP END PARALLEL DO

            !
            ! VH update
            !
            call linear_dispersive_cellcentred_matmult_JACOBI_ADAPTIVE(U(:,:,ELV), U(:,:,UH), U(:,:,VH), &
                msl_linear, distance_bottom_edge, distance_left_edge, dlon, dlat, &
                ds%offdiagonalAx(:,:,UH), ds%offdiagonalAx(:,:,VH), &
                ds%diagonalA(:,:,UH), ds%diagonalA(:,:,VH), &
                update_UH=.false., update_VH=.true., needs_update=ds%needs_update, i0=ds%i0, i1=ds%i1, &
                td1 = ds%td1, td2 = ds%td2)

            !$OMP PARALLEL DO SCHEDULE(DYNAMIC, 10) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            !!$OMP PARALLEL DO SCHEDULE(STATIC) DEFAULT(PRIVATE) SHARED(ds, U, msl_linear) REDUCTION(max: max_err)
            do j = 2, size(U,2) - 1
                !do i = 2, size(U, 1) - 1
                do i = ds%i0(j), ds%i1(j)
                    if(.not. ds%needs_update(i,j)) cycle

                    last_U = U(i,j,VH)
                    U(i,j,VH) = (ds%RHS(i,j,VH) + ds%offdiagonalAx(i,j,VH))/(1.0_dp - ds%diagonalA(i,j,VH))
                    U(i,j,VH) = U(i,j,VH) + jacobi_overrelax*(U(i,j,VH)-last_U)

                    ! Record the max abs_vh_difference/depth, reducing to abs_vh_difference in depths < 1m
                    ds%local_err(i,j) = max(ds%local_err(i,j), abs(U(i, j,VH) - last_U)/&
                         max(msl_linear - U(i,j,ELV), 1.0_dp))
                    max_err = max( max_err, ds%local_err(i,j) )
                end do
            end do
            !$OMP END PARALLEL DO

            ds%max_err = max_err

            ! Check for tolerance
            if(max_err < ds%tol) then
                exit jacobi_iter
                !if(allowed_to_exit) exit jacobi_iter
                !! Once the error is less than the tolerance, we 
                !! start doing full grid iterations again. 
                !! Ideally we still meet the tolerance after one, and ensure 
                !! the tolerance is met everywhere.
                !! 
                !allowed_to_exit = .true.
                !ds%i0 = 2
                !ds%i1 = size(U, 1)-1
                !!$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(ds)
                !do j = 1, size(U,2)
                !    ds%needs_update(:,j) = .true.
                !end do
                !!$OMP END PARALLEL DO
            end if

            if(iter > (min_jacobi_iterations - 1) .and. (.not. allowed_to_exit)) then
                ! Determine which areas need updates.
                ! We record cells if cells are 'near' a cell that has error > "some fraction of tol" (ds%needs_update(:,:)).
                ! We also record, for each j, a range of i-indices that include all those needing updates.
                !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(ds)

                !$OMP DO SCHEDULE(DYNAMIC, 10)
                do j = 2, size(U,2) - 1
                    ! Update the 'needs_update' variable

                    ! We don't need to scan every i-index in the j'th row
                    ! Need to look look at range of active cells 'bw' rows above/below the current row,
                    ! as these might have error > tol. Need care to not go out of bounds. 
                    l0 = max(2, j-bw) ! safe j index below
                    l1 = min(size(U,2)-1, j+bw) ! safe j index above

                    ! i indices that we should look at in the j'th row
                    ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
                    ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above

                    do i = ii0 , ii1 
                        ! Check within local window for error > tol
                        i0 = max(1, i-bw)
                        i1 = min(size(U,1), i+bw)
                        j0 = max(1, j-bw)
                        j1 = min(size(U,2), j+bw)
                        ds%needs_update(i, j) = any(ds%local_err(i0:i1,j0:j1) > (ds%tol*ds%local_tol_factor))
                    end do
                end do
                !$OMP END DO

                !$OMP DO SCHEDULE(DYNAMIC, 10)
                do j = 2, size(U, 2) - 1
                    ! Update the 'ds%i0, ds%i1' variables defining the i-range in which
                    ! we search for each j. Could not do this in the previous loop due to
                    ! race condition.

                    ! See the loop above for documentation of the next 4 lines
                    l0 = max(2, j-bw) ! safe j index below
                    l1 = min(size(U,2)-1, j+bw) ! safe j index above
                    ii0 = max(2, minval(ds%i0(l0:l1) - bw)) ! Safe i-index below
                    ii1 = min(size(U,1)-1, maxval(ds%i1(l0:l1) + bw)) ! Safe i-index above

                    ! To update in loop -- 
                    ! i0 = min index to update; initialise with a large value
                    ds%i0(j) = size(U, 1)
                    ! i1 = max index to update; initialise with a small value
                    ds%i1(j) = 1
                    do i = ii0, ii1
                        if(ds%needs_update(i,j)) then
                            ! Cell i,j is at least close to a cell that has error above tolerance.
                            ds%i0(j) = min(ds%i0(j), i)
                            ds%i1(j) = max(ds%i1(j), i)
                        end if
                    end do
                end do
                !$OMP END PARALLEL
            end if

        end do jacobi_iter

        !print*, 'Jacobi iter: ', ds%last_iter, max_err, ', shape(U): ', shape(U) !, count(ds%needs_update)*1.0_dp/(size(U,1)*size(U,2))
        !if(ds%last_iter == ds%max_iter) then
        !    write(log_output_unit, *) 'Jacobi iteration hit max iterations (', ds%max_iter, ') with error ', max_err
        !end if

    end subroutine
#endif

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

                ! See comments at the start of this code for explanation of these formulas, which apply in both spherical and
                ! cartesian coordinates.
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

    subroutine linear_dispersive_solve_staggered_grid_TRIDIAG(ds, U, &
            dlon, dlat, distance_bottom_edge, distance_left_edge, msl_linear)
        class(dispersive_solver_type), intent(inout) :: ds
        real(dp), intent(inout) :: U(:,:,:)
            ! domain%U after an explicit shallow water update
        real(dp), intent(in) :: dlon, dlat, distance_bottom_edge(:), distance_left_edge(:), msl_linear
            ! cell dx, dy, edge distances, and mean-sea-level for the linear solver

        integer :: i, j, iter
        real(dp) :: max_err, last_U

        ! Get the explicit part of the dispersive terms
        call linear_dispersive_staggered_matmult_JACOBI(&
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
