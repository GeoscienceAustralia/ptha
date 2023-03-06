module local_routines 
    !!
    !! Setup the Merewether problem.
    !!

    use global_mod, only: dp, ip, force_double, charlen, pi
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: multi_raster_type
    use which_mod, only: which
    use ragged_array_mod, only: ragged_array_2d_ip_type
    use file_io_mod, only: read_character_file, read_csv_into_array
    use points_in_poly_mod, only: points_in_poly
    use burn_into_grid_mod, only: polygons_values_type
    use logging_mod, only: log_output_unit
    use iso_c_binding, only: c_f_pointer, c_loc
    implicit none

    ! Add rain
    ! Discharge of 19.7 m^3/s occurs over a cirle with radius 15
    real(dp), parameter :: rain_centre(2) = [382300.0_dp, 6354290.0_dp], rain_radius = 15.0_dp
    real(force_double) :: rain_rate =  real(19.7_dp, force_double)/(pi * rain_radius**2)
    integer(ip) :: num_input_discharge_cells

    contains 

    !
    ! Main setup routine
    !
    subroutine set_initial_conditions_merewether(domain)            
        class(domain_type), target, intent(inout):: domain

        integer(ip):: i, j, k, house_file_unit
        real(dp), allocatable:: x(:), y(:)

        ! For reading background elevation
        character(len=charlen):: input_elevation_file
        type(multi_raster_type) :: elevation_data

        ! Add this elevation to all points in houses/ csv polygons
        real(dp), parameter:: house_height = 3.0_dp
        ! friction parameters
        real(dp), parameter:: friction_road = 0.02_dp, friction_other = 0.04_dp

        ! For setting values in polygons
        character(len=charlen), allocatable:: house_filenames(:)
        real(dp), allocatable :: house_heights(:)
        type(polygons_values_type) :: house_pvt, road_pvt
        
        ! For inflow boundary
        type(ragged_array_2d_ip_type), pointer :: rainfall_region_indices


        ! Dry flow to start with. Later we clip stage to elevation.
        domain%U(:,:,[STG, UH, VH]) = 0.0_dp

        ! Set elevation with the raster
        input_elevation_file = "./topography/topography1.tif"
        call elevation_data%initialise([input_elevation_file])

        ! To hold x/y coordinates along a row with constant y-value
        allocate(x(domain%nx(1)), y(domain%nx(1)))
        x = domain%x

        ! Read the elevation data row-by-row
        do j = 1, domain%nx(2)
            y = domain%y(j)
            call elevation_data%get_xy(x, y, domain%U(:,j,ELV), domain%nx(1), bilinear=0_ip)
        end do
        call elevation_data%finalise

        write(log_output_unit, *) 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        ! Add house heights to DEM
        open(newunit = house_file_unit, file='houses_filenames.txt') ! Names of house polygons
        call read_character_file(house_file_unit, house_filenames, "(A)")
        close(house_file_unit)
        allocate(house_heights(size(house_filenames)))
        house_heights = house_height ! Height in each house
        call house_pvt%setup(house_filenames, house_heights) ! Type for setting grid values in polygons
        call house_pvt%burn_into_grid(domain%U(:,:,ELV), domain%lower_left, domain%lower_left+domain%lw, burn_type='add')

        ! Transmissive outflow
        ! We prevent mass leaking by using walls for 2 boundaries
        domain%U(:, 1:2, ELV) = 100.0_dp
        domain%U(1:2, :, ELV) = 100.0_dp
        ! We use a flat boundary on the east and north so that the transmissive BC is well-behaved
        domain%U(domain%nx(1), :, ELV) = domain%U(domain%nx(1)-1, :, ELV)
        domain%U(:, domain%nx(2), ELV) = domain%U(:, domain%nx(2)-1, ELV) 

        ! Set friction
        domain%manning_squared = friction_other * friction_other ! Background value
        ! Another friction value for points on the road
        call road_pvt%setup([character(len=charlen)::'Road/RoadPolygon.csv'], [friction_road*friction_road])
        call road_pvt%burn_into_grid(domain%manning_squared, domain%lower_left, domain%lower_left+domain%lw)

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-08_dp)

        write(log_output_unit, *) 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
        write(log_output_unit, *) 'Stage range: ', minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))

        ! Figure out indices that are inside the rainfall forcing region.
        allocate(rainfall_region_indices)
        allocate(rainfall_region_indices%i2(domain%nx(2)))
        num_input_discharge_cells = 0
        do j = 1, domain%nx(2)
            ! Find cells in this row that are within the rain circle
            call which( (domain%x - rain_centre(1))**2 < rain_radius**2 - (domain%y(j) - rain_centre(2))**2, &
                rainfall_region_indices%i2(j)%i1)

            ! For this routine we cross-check conservation by counting the number of input discharge cells and
            ! doing a separate mass balance. This assumes 1 domain (only)
            num_input_discharge_cells = num_input_discharge_cells + size(rainfall_region_indices%i2(j)%i1)
        end do

        ! Add rainfall to the domain
        domain%forcing_subroutine => apply_rainfall_forcing
        domain%forcing_context_cptr = c_loc(rainfall_region_indices)
        call domain%store_forcing()

        deallocate(x,y, house_filenames)

    end subroutine

    !
    ! This is the discharge source term, called every time-step. It's like rainfall 
    ! in a "circle"
    !
    subroutine apply_rainfall_forcing(domain, dt)
        type(domain_type), intent(inout) :: domain
        real(dp), intent(in) :: dt

        integer(ip) :: j
        type(ragged_array_2d_ip_type), pointer :: rainfall_region_indices

        ! Unpack the forcing context pointer. In this case it is a ragged array giving the indices
        ! where we should apply the rainfall in this domain
        call c_f_pointer(domain%forcing_context_cptr, rainfall_region_indices)

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(domain, rain_rate, dt, rainfall_region_indices)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            if(size(rainfall_region_indices%i2(j)%i1) > 0) then
                ! Rainfall at cells within a circle of radius "rain_radius" about "rain_centre"
                domain%U(    rainfall_region_indices%i2(j)%i1, j,STG) = rain_rate * dt + &
                    domain%U(rainfall_region_indices%i2(j)%i1, j,STG) 
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine


end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program merewether
    !! Urban flooding test case in Merewether from Australian Rainfall and Runoff.
    !! Smith, G. & Wasko, C. Revision Project 15: Two Dimensional Simulations In 
    !! Urban Areas - Representation of Buildings in 2D Numerical Flood Models Australian 
    !! Rainfall and Runoff, Engineers Australia, 2012

    use global_mod, only: ip, dp, minimum_allowed_depth, default_nonlinear_timestepping_method
    use multidomain_mod, only: multidomain_type
    use boundary_mod, only: transmissive_boundary
    use local_routines

    implicit none

    integer(ip):: j
    real(dp):: last_write_time
    type(multidomain_type):: md

    ! Global timestep
    real(dp), parameter :: global_dt = 0.08_dp

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 60.0_dp
    real(dp), parameter :: final_time = 900.0_dp
    
    ! Length/width
    real(dp), parameter, dimension(2):: global_lw = [320.0_dp, 415.0_dp]
    ! Lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = [382251.0_dp, 6354266.5_dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [320_ip, 415_ip] ![160_ip, 208_ip] ![321_ip, 416_ip]
    ! Track mass
    real(force_double) :: volume_added

    ! Single domain model 
    allocate(md%domains(1)) 

    md%domains(1)%lw = global_lw
    md%domains(1)%nx = global_nx
    md%domains(1)%lower_left = global_ll
    md%domains(1)%timestepping_method = default_nonlinear_timestepping_method

    ! Initialise the domain
    call md%setup()

    do j = 1, size(md%domains)
        ! Set the initial and boundary conditions
        call set_initial_conditions_merewether(md%domains(j))

        if(.not. all(md%domains(j)%is_nesting_boundary)) then
            ! Allow waves to propagate out of the domain 
            md%domains(j)%boundary_subroutine => transmissive_boundary
            ! Note we adjusted the topography 
            ! (inside set_initial_conditions_merewether) to ensure that
            ! the outflow is restricted to the downstream edge, and that the
            ! transmissive boundary is well behaved (boundary elevation 
            ! equal to that just inside the boundary).
        end if
    end do

    call md%make_initial_conditions_consistent() 

    ! Evolve the model
    do while (.TRUE.)

        call md%write_outputs_and_print_statistics(approximate_writeout_frequency=approximate_writeout_frequency)

        if (md%domains(1)%time > final_time) exit

        call md%evolve_one_step(global_dt)

    end do

    write(log_output_unit, *) ''
    write(log_output_unit, '(A)') 'Expected mass change due to inflows (that SWALS does not account for in mass-tracking):'
    write(log_output_unit, *) num_input_discharge_cells * md%domains(1)%time * product(md%domains(1)%dx) * rain_rate
    write(log_output_unit, *) ''

    call md%finalise_and_print_timers()


END PROGRAM
