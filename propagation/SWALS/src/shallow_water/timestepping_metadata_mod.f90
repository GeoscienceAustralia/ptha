module timestepping_metadata_mod
    !!
    !! Keep most solver metadata here in an array of type "timestepping_metadata_type", 
    !! where each entry contains important metadata for one numerical method.
    !!

    use global_mod, only : dp, ip, charlen
    use logging_mod, only : log_output_unit
    use stop_mod

    implicit none

    ! Number of timestepping methods
    integer, parameter, private :: n_ts = 8

    type timestepping_metadata_type
        !! Type that holds metadata for each type of solver
        character(len=charlen) :: timestepping_method = '' 
            !! Name of the timestepping method
        
        integer(ip) :: is_staggered_grid = 0 
            !! Is the grid treated as staggered (1) or colocated (0) ?
        logical :: flux_correction_is_unsupported = .false. 
            !! Flag for solvers that do not track fluxes (so cannot do any flux-correction)
        logical :: flux_correction_of_mass_only = .false. 
            !! Some solvers can only do flux correction of mass, but not momentum
        real(dp) :: default_cfl = -1.0_dp
            !! CFL condition
        real(dp) :: default_theta = -1.0_dp
            !! Parameter affecting the slope-limiter for finite-volume methods.
        logical :: adaptive_timestepping = .true.
            !! Can the solver compute its own time-step adaptively?
       
        integer(ip) :: nesting_thickness_for_one_timestep = -1_ip
            !! How many halo cells are required to advance interior cells a single timestep while
            !! retaining a valid solution?
            !! For example, consider the 'euler' finite-volume timestepping method.  A single 'euler' step of the finite volume solver
            !! needs a 2-layer halo. Suppose the cells are indexed 1, 2, 3, ... N, and initially contain valid values. If we try to
            !! evolve in time by one step, then 'cell 1' has no-way to compute the flux at 'edge (1-1/2)' -- and furthermore, the flux at
            !! 'edge (1+1/2)' is problematic because 'cell 1' cannot compute its gradient, so the update of 'cell 2' is also invalid. 
            !! Thus, a 2-layer halo is needed (where the values for 'cell 1', 'cell 2', 'cell (N-1)', 'cell N' are provided
            !! by halo exchanges with another grid, or by boundary-condition assumptions).

        character(len=charlen) :: forcing_elevation_is_allowed = 'never'
            !! SWALS supports user-provided forcing terms (via domain%forcing_subroutine). This can in-principle also
            !! operate on the domain elevation. However, some solvers cannot do this (or need special options enabled to
            !! do it). Allowed values are "always" [for schemes that don't have any issues] and "optional" [for schemes
            !! where it must be enabled manually, by setting domain%elevation_forcing_allowed = .TRUE.], and "never".
    end type

    type(timestepping_metadata_type), public, protected :: timestepping_metadata(n_ts)
        !!
        !! Main source of metadata. Initialise this in a subroutine (although it's static - ideally this would
        !! be a parameter -- not sure how to do that neatly in Fortran)
        !!

    logical, private :: IS_SETUP = .FALSE.
        !! Record whether we've setup the metadata, so we don't have to call setup_timestepping_metadata too often.

    contains

        subroutine setup_timestepping_metadata()
            !! Populate timestepping_metadata.
            !! Note this only needs to be called once

            !
            ! rk2 defaults
            !
            timestepping_metadata(1)%timestepping_method = 'rk2'
            ! NOTE: In some case rk2 is stable with a larger time-step, see 
            ! Andrew Giuliani, Lilia Krivodonova, On the optimal CFL
            !   number of SSP methods for hyperbolic problems.
            timestepping_metadata(1)%default_cfl = 0.99_dp
            timestepping_metadata(1)%default_theta = 1.6_dp
            timestepping_metadata(1)%nesting_thickness_for_one_timestep = 4_ip
            timestepping_metadata(1)%forcing_elevation_is_allowed = 'optional'

            !
            ! rk2n defaults
            !
            timestepping_metadata(2)%timestepping_method = 'rk2n'
            timestepping_metadata(2)%default_cfl = 0.99_dp
            timestepping_metadata(2)%default_theta = 1.6_dp
            timestepping_metadata(2)%nesting_thickness_for_one_timestep = 10_ip
            timestepping_metadata(2)%forcing_elevation_is_allowed = 'optional'

            !
            ! midpoint defaults
            !
            timestepping_metadata(3)%timestepping_method = 'midpoint'
            timestepping_metadata(3)%default_cfl = 0.99_dp
            timestepping_metadata(3)%default_theta = 1.6_dp
            timestepping_metadata(3)%nesting_thickness_for_one_timestep = 4_ip
            timestepping_metadata(3)%forcing_elevation_is_allowed = 'optional'

            !
            ! Euler defaults
            !
            timestepping_metadata(4)%timestepping_method = 'euler'
            timestepping_metadata(4)%default_cfl = 0.9_dp
            timestepping_metadata(4)%default_theta = 0.9_dp
            timestepping_metadata(4)%nesting_thickness_for_one_timestep = 2_ip
            timestepping_metadata(4)%forcing_elevation_is_allowed = 'always'

            !
            ! Linear leapfrog defaults
            !
            timestepping_metadata(5)%timestepping_method = 'linear'
            timestepping_metadata(5)%is_staggered_grid = 1
            ! Updating the momentum advection flux between linear/nonlinear can cause instabilities, unsurprisingly
            ! given that this flux is not included in the linear equations
            timestepping_metadata(5)%flux_correction_of_mass_only = .true.
            timestepping_metadata(5)%default_cfl = 0.7_dp
            timestepping_metadata(5)%nesting_thickness_for_one_timestep = 2_ip
            timestepping_metadata(5)%adaptive_timestepping = .false.
            timestepping_metadata(5)%forcing_elevation_is_allowed = 'always'

            !
            ! Linear leapfrog + nonlinear friction defaults
            !
            timestepping_metadata(6)%timestepping_method = 'leapfrog_linear_plus_nonlinear_friction'
            timestepping_metadata(6)%is_staggered_grid = 1
            ! Updating the momentum advection flux between linear/nonlinear can cause instabilities, unsurprisingly
            ! given that this flux is not included in the linear equations
            timestepping_metadata(6)%flux_correction_of_mass_only = .true.
            timestepping_metadata(6)%default_cfl = 0.7_dp
            timestepping_metadata(6)%nesting_thickness_for_one_timestep = 2_ip
            timestepping_metadata(6)%adaptive_timestepping = .false.
            timestepping_metadata(6)%forcing_elevation_is_allowed = 'always'

            !
            ! Cliffs (by Elena Tolkova -- similar to MOST)
            !
            timestepping_metadata(7)%timestepping_method = 'cliffs'
            timestepping_metadata(7)%is_staggered_grid = 0
            ! No flux correction for cliffs -- because this is not yet implemented, and in any case the
            ! solver is not mass conservative
            timestepping_metadata(7)%flux_correction_is_unsupported = .true.
            timestepping_metadata(7)%default_cfl = 0.7_dp
            timestepping_metadata(7)%nesting_thickness_for_one_timestep = 2_ip
            timestepping_metadata(7)%forcing_elevation_is_allowed = 'always'


            !
            ! Nonlinear leapfrog with upwind advection
            !
            timestepping_metadata(8)%timestepping_method = 'leapfrog_nonlinear'
            timestepping_metadata(8)%is_staggered_grid = 1
            !timestepping_metadata(8)%flux_correction_of_mass_only = .true.
            timestepping_metadata(8)%default_cfl = 0.7_dp
            timestepping_metadata(8)%nesting_thickness_for_one_timestep = 3_ip
            timestepping_metadata(8)%adaptive_timestepping = .false.
            timestepping_metadata(8)%forcing_elevation_is_allowed = 'always'

        end subroutine

        function timestepping_method_index(timestepping_method) result(ts_index)
            !! Given a timestepping_method (e.g. 'linear' or 'rk2'),
            !! return the corresponding index of timestepping_metadata
            character(len=*), intent(in) :: timestepping_method 
                !! The timestepping method (e.g. 'rk2' or 'linear', etc)
            integer(ip) :: ts_index

            integer(ip) :: i

            if(.not. IS_SETUP) then
                call setup_timestepping_metadata
                IS_SETUP = .TRUE.
            end if

            ts_index = -HUGE(1_ip) ! Unset

            ! Find the index of timestepping_metadata matching the timestepping method
            do i = 1, n_ts
                if(TRIM(timestepping_metadata(i)%timestepping_method) == TRIM(timestepping_method)) ts_index = i
            end do

            if(ts_index == -HUGE(1_ip)) then
                write(log_output_unit,*) 'timestepping_method: ', TRIM(timestepping_method), ' not recognized'
                call generic_stop
            end if

        end function

end module
