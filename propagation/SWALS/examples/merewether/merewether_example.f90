module local_routines 
    !!
    !! Setup the Merewether problem.
    !!

    use global_mod, only: dp, ip, charlen, wall_elevation
    use domain_mod, only: domain_type, STG, UH, VH, ELV
    use read_raster_mod, only: read_gdal_raster
    use which_mod, only: which
    use ragged_array_mod, only: ragged_array_2d_ip_type
    use file_io_mod, only: read_character_file, read_csv_into_array, count_file_lines
    use points_in_poly_mod, only: points_in_poly
    implicit nONE

    contains 

    !
    ! Main setup routine
    !
    subroutine set_initial_conditions_merewether(domain, input_discharge_indices, reflective_boundaries)            

        class(domain_type), target, intent(inout):: domain
        type(ragged_array_2d_ip_type), intent(inout) :: input_discharge_indices
        logical, intent(in):: reflective_boundaries

        integer(ip):: i, j, k
        real(dp), allocatable:: x(:,:), y(:,:)

        character(len=charlen):: input_elevation_file, polygon_filename

        ! Add this elevation to all points in houses/ csv polygons
        real(dp), parameter:: house_height = 3.0_dp
        ! friction parameters
        real(dp), parameter:: friction_road = 0.02_dp, friction_other = 0.04_dp

        ! things to help read polygons
        character(len=charlen), allocatable:: house_filenames(:)
        integer(ip):: house_file_unit, inside_point_counter, house_cell_count
        real(dp), allocatable:: polygon_coords(:,:)
        logical, allocatable:: is_inside_poly(:)


        allocate(x(domain%nx(1),domain%nx(2)), y(domain%nx(1),domain%nx(2)))

        !
        ! Stage
        !

        domain%U(:,:,[STG, UH, VH]) = 0.0_dp

        !
        ! Set elevation with the raster
        !

        input_elevation_file = "./topography/topography1.tif"

        do j = 1, domain%nx(2)
          do i = 1, domain%nx(1)
            x(i,j) = domain%lower_left(1) + (i-0.5_dp)*domain%dx(1) 
            y(i,j) = domain%lower_left(2) + (j-0.5_dp)*domain%dx(2) 
          end dO
        end do

        call read_gdal_raster(input_elevation_file, x, y, domain%U(:,:,ELV), &
            domain%nx(1)*domain%nx(2), verbose=1_ip, bilinear=0_ip)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))

        ! Get filenames for the houses
        open(newunit = house_file_unit, file='houses_filenames.txt')
        call read_character_file(house_file_unit, house_filenames, "(A)")
        close(house_file_unit)

        ! Read the houses and add house_height to the elevation for all points inside
        allocate(is_inside_poly(domain%nx(1)))
        house_cell_count = 0
        do k = 1, size(house_filenames)

            inside_point_counter = 0
            polygon_filename = './' // TRIM(house_filenames(k))
            call read_csv_into_array(polygon_coords, polygon_filename)

            do j = 1, domain%nx(2)
                ! Find points in the polygon in column j
                call points_in_poly(polygon_coords(1,:), polygon_coords(2,:), x(:,j), y(:,j), is_inside_poly)
                inside_point_counter = inside_point_counter + count(is_inside_poly)

                do i = 1, domain%nx(1) 
                    if(is_inside_poly(i)) then
                        domain%U(i,j,ELV) = domain%U(i,j,ELV) + house_height
                        house_cell_count = house_cell_count + 1_ip
                    end if
                end do

            end do 

            print*, trim(house_filenames(k)), ': ', inside_point_counter

        end do

        print*, '# cells in houses: ', house_cell_count
       
        if(reflective_boundaries) then 
            domain%U(1,:,ELV) = wall_elevation
            domain%U(domain%nx(1), :, ELV) = wall_elevation    
            domain%U(:, 1 ,ELV) = wall_elevation
            domain%U(:, domain%nx(2), ELV) = wall_elevation    
        END IF

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-06_dp)
        
        !
        ! Set friction
        !

        ! Initial value
        domain%manning_squared = friction_other * friction_other

        ! Find points in road polygon
        polygon_filename = 'Road/RoadPolygon.csv'
        call read_csv_into_array(polygon_coords, polygon_filename)
        inside_point_counter = 0
        do j = 1, domain%nx(2)

            call points_in_poly(polygon_coords(1,:), polygon_coords(2,:), x(:,j), y(:,j), is_inside_poly)
            inside_point_counter = inside_point_counter + count(is_inside_poly)

            do i = 1, domain%nx(1) 
                if(is_inside_poly(i)) domain%manning_squared(i,j) = friction_road * friction_road
            end do

        end do

        print*, ''
        print*, '# Points in road polygon :', inside_point_counter
        print*, ''

        !       
        ! Get input discharge indices 
        !
        ! The anuga code was:
        ! fixed_inflow = anuga.Inflow(domain,
        !           center=(382300.0,6354290.0),
        !           radius=15.00,
        !           rate=19.7)

        allocate(input_discharge_indices%i2(domain%nx(2)))

        ! For each column, find the indices within 15m of a particular point
        do j = 1, domain%nx(2)
            call which( (x(:,j) - 382300.0)**2 + (y(:,j) - 6354290.0)**2 < 15.0_dp**2, &
                input_discharge_indices%i2(j)%i1)
        end do

        if(.not. reflective_boundaries) THEN
            ! prevent mass leaking out of the back of the domain
            domain%U(:, 1:2, ELV) = domain%U(:,1:2, ELV) + 10.0_dp
        end if 

        ! Ensure stage >= elevation
        domain%U(:,:,STG) = max(domain%U(:,:,STG), domain%U(:,:,ELV) + 1.0e-06_dp)

        print*, 'Elevation range: ', minval(domain%U(:,:,ELV)), maxval(domain%U(:,:,ELV))
        print*, 'Stage range: ', minval(domain%U(:,:,STG)), maxval(domain%U(:,:,STG))

    end subroutine

end module 

!@!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program merewether
    !! Urban flooding test case in Merewether from Australian Rainfall and Runoff.
    !! Smith, G. & Wasko, C. Revision Project 15: Two Dimensional Simulations In 
    !! Urban Areas - Representation of Buildings in 2D Numerical Flood Models Australian 
    !! Rainfall and Runoff, Engineers Australia, 2012

    use global_mod, only: ip, dp, minimum_allowed_depth
    use ragged_array_mod, only: ragged_array_2d_ip_type
    use domain_mod, only: domain_type
    use boundary_mod, only: transmissive_boundary
    use local_routines
    use iso_c_binding, only: C_DOUBLE

    implicit none

    integer(ip):: j
    real(dp):: last_write_time, rain_rate, mass_integral, domain_volume
    real(C_DOUBLE) :: volume_added
    type(domain_type):: domain
    real(dp), parameter:: pi = atan(1.0_dp) * 4.0_dp

    ! Approx timestep between outputs
    real(dp), parameter :: approximate_writeout_frequency = 60.0_dp
    real(dp), parameter :: final_time = 900.0_dp
    
    ! Length/width
    real(dp), parameter, dimension(2):: global_lw = [320.0_dp, 415.0_dp]
    ! Lower-left corner coordinate
    real(dp), parameter, dimension(2):: global_ll = [382251.0_dp, 6354266.5_dp]
    ! grid size (number of x/y cells)
    integer(ip), parameter, dimension(2):: global_nx = [320_ip, 415_ip] ![160_ip, 208_ip] ![321_ip, 416_ip]
    ! Use reflective or transmissive boundaries
    logical, parameter:: reflective_boundaries = .FALSE.
   
    ! indices where the input discharge goes 
    type(ragged_array_2d_ip_type):: input_discharge_indices
    integer(ip) :: num_input_discharge_indices
    
    !domain%theta = 1.00_dp
    domain%timestepping_method = 'rk2n' !'euler' !'rk2n'
    domain%maximum_timestep = 5.0_dp
    
    ! Allow waves to propagate outside all edges
    if(.not. reflective_boundaries) then
        domain%boundary_subroutine => transmissive_boundary
    end if

    ! Allocate domain -- must have set timestepping method BEFORE this
    call domain%allocate_quantities(global_lw, global_nx, global_ll)

    ! Call local routine to set initial conditions
    call set_initial_conditions_merewether(domain, input_discharge_indices, reflective_boundaries)
    call domain%update_boundary() ! Make boundary consistent
    volume_added = domain%volume_interior()
    
    ! Trick to get the code to write out just after the first timestep
    last_write_time = -approximate_writeout_frequency

    ! Count the input discharge indices
    num_input_discharge_indices = 0
    do j = 1, size(input_discharge_indices%i2)
        num_input_discharge_indices = num_input_discharge_indices + size(input_discharge_indices%i2(j)%i1)
        !print*, ' ## ', j, ' ## ', input_discharge_indices%i2(j)%i1
    end do
    print*, '# discharge indices :', num_input_discharge_indices

    ! Evolve the code
    do while (.TRUE.)

        if(domain%time - last_write_time >= approximate_writeout_frequency) then

            last_write_time = last_write_time + approximate_writeout_frequency

            call domain%print()
            call domain%write_to_output_files()

            print*, 'mass_balance: ', domain%mass_balance_interior() 
            print*, 'volume_added: ', volume_added
            print*, '.... difference: ', domain%mass_balance_interior() - volume_added
        end if

        if (domain%time > final_time) then
            exit 
        end if

        call domain%evolve_one_step()

        !print*, 'MAX-dt: ', domain%max_dt

        ! Add rain
        ! Discharge of 19.7 m^3/s occurs over a cirle with radius 15
        rain_rate = domain%evolve_step_dt * (19.7_dp/(pi * 15.0_dp**2))

        volume_added = volume_added + rain_rate * num_input_discharge_indices * product(domain%dx)

        !$OMP PARALLEL SHARED(domain, rain_rate, input_discharge_indices)
        !$OMP DO SCHEDULE(STATIC)
        do j = 1, domain%nx(2)
            if (size(input_discharge_indices%i2(j)%i1) > 0) then
                domain%U(input_discharge_indices%i2(j)%i1,j,1) = &
                domain%U(input_discharge_indices%i2(j)%i1,j,1) + rain_rate
            end if
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end do

    call domain%write_max_quantities()

    ! Print timing info
    call domain%timer%print()

    ! Simple test (mass conservation only)
    print*, 'MASS CONSERVATION TEST'
    if(abs(domain%mass_balance_interior() - volume_added) < 1.0e-06_dp) then
        print*, 'PASS'
    else
        print*, 'FAIL -- mass conservation not good enough. ', &
            'This is expected with single precision, which will affect the computed direct rainfall volumes'
    end if

    call domain%finalise()
END PROGRAM
