MODULE point_gauge_mod
    !
    ! Type to store and write out gauge time-series at chosen locations
    !
    USE global_mod, only: dp, ip, output_precision, charlen
    USE file_io_mod, only: read_csv_into_array
    USE stop_mod, only: generic_stop

    ! point_gauge output can be either:
    !   A) netcdf, or 
    !   B) a home-brew mix of ascii and binary files
    ! The latter can be turned on by compiling with -DNONETCDF. 
    ! This is mainly useful if it is difficult to link with netcdf (e.g. due to a different
    ! compiler having been used to compile netcdf vs the program)
#ifndef NONETCDF
    USE netcdf
#endif

    IMPLICIT NONE
    
    PRIVATE
    PUBLIC:: point_gauge_type, test_point_gauge_mod

    ! Note -- these definitions must agree with those in domain_mod.f90
    ! We don't import them here since domain_mod already imports the point_gauge_mod
    INTEGER(ip), PARAMETER :: STG=1, UH=2, VH=3, ELV=4

#ifndef NONETCDF
    ! If we don't close the netcdf file, we have to flush
    ! it to ensure all writes occur. Probably inefficient
    !LOGICAL, PARAMETER :: flush_netcdf=.FALSE.
#endif

    ! Used to store values from an array at an irregular set of points
    !
    ! The array from which point gauges are extracted is assumed to be 3 dimensional U(:,:,:)
    ! The first 2 dimensions are spatial (x/y or lon/lat), and 
    ! the last dimension represents the quantity (e.g. Stage, UH, VH, ELEV)
    TYPE point_gauge_type

#ifdef NONETCDF
        !! Variables for writing binary/ascii output
        CHARACTER(charlen):: static_output_file, time_series_output_file, gauge_metadata_file
        INTEGER(ip):: static_output_unit, time_series_output_unit, gauge_metadata_unit
#else
        !! Variables for writing out to netcdf
        CHARACTER(charlen):: netcdf_gauge_output_file
        INTEGER :: netcdf_gauge_output_file_ID 
        INTEGER(ip) :: netcdf_time_var_ID
        INTEGER(ip), ALLOCATABLE:: time_series_ncdf_iVar_ID(:)
        ! Count how many steps have been output
        INTEGER(ip) :: netcdf_num_output_steps = 0
#endif

        ! How many gauges, and the dimension of the space
        INTEGER(ip):: n_gauges, space_dim = 2
        ! Indices of variables 
        INTEGER(ip), ALLOCATABLE:: time_series_var(:), static_var(:)

        ! Coordinates of gauges
        REAL(dp), ALLOCATABLE:: xy(:,:)
        ! Integer ID of gauges
        REAL(dp), ALLOCATABLE:: gauge_ids(:)
        ! Indices used to lookup gauge locations in domain
        INTEGER(ip), ALLOCATABLE:: site_index(:,:)
        ! Values stored at gauges each time
        REAL(dp), ALLOCATABLE:: time_series_values(:,:)
        ! Values stored at gauges when it is created
        REAL(dp), ALLOCATABLE:: static_values(:,:)

        CONTAINS

        PROCEDURE :: allocate_gauges => allocate_gauges
        PROCEDURE :: initialise_gauges => initialise_gauges
        PROCEDURE :: update_gauge_var => update_gauge_var
        PROCEDURE :: write_current_time_series => write_current_time_series
        PROCEDURE :: finalise => finalise_point_gauges

    END TYPE

    CONTAINS

    !
    ! Read gauge xy and allocate space for their data
    !
    ! @param point_gauges point_gauge_type
    ! @param xy_coordinates array of rank 2 with size [2, n_gauges] giving x,y
    !  coordinates of gauges
    ! @param time_series_var_indices Indices of variables in domain%U that we
    !  wish to store each time-step. e.g. [STG, UH, VH], or [STG]
    ! @param static_var_indices Indices of variables in domain%U that we only
    !  wish to store once. e.g. [ELV], or [UH, VH, ELV]
    ! @param gauge_ids real array of size n_gauges giving an ID for each
    !  gauge. Make it REAL to avoid truncation issues for large IDs
    !
    SUBROUTINE allocate_gauges(point_gauges, xy_coordinates, &
        time_series_var_indices, static_var_indices, gauge_ids)
        CLASS(point_gauge_type), INTENT(INOUT) :: point_gauges
        REAL(dp), INTENT(IN) :: xy_coordinates(:,:)
        INTEGER(ip), INTENT(IN) :: time_series_var_indices(:), static_var_indices(:)
        REAL(dp), INTENT(IN) :: gauge_ids(:)

        INTEGER(ip):: n_gauges, space_dim, n_ts_var, n_static_var
            
            ! At the moment the space dimension must be 2 
            if (point_gauges%space_dim /= size(xy_coordinates(:,1))) then
                print*, 'gauge xy dimension is not equal ', &
                    point_gauges%space_dim
            end if

            ! Get coordinates
            n_gauges = size(xy_coordinates(1,:))
            space_dim = point_gauges%space_dim
            point_gauges%n_gauges = n_gauges
            allocate(point_gauges%xy(space_dim, n_gauges))
            point_gauges%xy = xy_coordinates 
            
            space_dim = point_gauges%space_dim
            n_gauges = point_gauges%n_gauges

            n_ts_var = size(time_series_var_indices)
            n_static_var = size(static_var_indices)
   
            ALLOCATE(point_gauges%site_index(space_dim, n_gauges))
            ALLOCATE(point_gauges%time_series_values(n_gauges, n_ts_var))
            ALLOCATE(point_gauges%static_values(n_gauges, n_static_var))
            ALLOCATE(point_gauges%gauge_ids(n_gauges))
            point_gauges%gauge_ids = gauge_ids

            ALLOCATE(point_gauges%time_series_var(n_ts_var))
            point_gauges%time_series_var = time_series_var_indices
            ALLOCATE(point_gauges%static_var(n_static_var))
            point_gauges%static_var = static_var_indices

    END SUBROUTINE

    ! Convenience subroutine to lookup variables in domain_U at point_gauges,
    ! and populate time_series_values with the values of quantities var_inds
    !
    ! This can be used to populate either point_gauges%static_values or 
    ! point_gauges%time_series_values
    !
    ! @param point_gauges point_gauge_type
    ! @param time_series_values array with size [n_gauges, size(var_inds)] that holds the time series data
    ! @param var_inds Array giving indices in the 3rd dimension of domain_U that we write out. e.g. [STG, UH, VH], or [STG]
    !
    SUBROUTINE update_gauge_var(point_gauges, time_series_values, domain_U, &
        var_inds)
        CLASS(point_gauge_type), INTENT(INOUT):: point_gauges
        REAL(dp), INTENT(INOUT):: time_series_values(:,:)
        REAL(dp), INTENT(IN):: domain_U(:,:,:)
        INTEGER(ip), INTENT(IN):: var_inds(:)

        INTEGER(ip):: i, j, xind, yind

        if(size(var_inds) /= size(time_series_values(1,:))) then
            print*, 'Dimension mismatch between time_series_values and var_inds'
            call generic_stop()
        end if

        DO j = 1, size(var_inds)
            DO i = 1, point_gauges%n_gauges
                xind = point_gauges%site_index(1,i)
                yind = point_gauges%site_index(2,i)
                time_series_values(i,j) = domain_U(xind, yind, var_inds(j))
            END DO 
        END DO

    END SUBROUTINE
      
    !
    ! Write time snapshot at gauges to netcdf
    ! 
    ! @param point_gauges point_gauges_type
    ! @param domain_U array from which we extract data at the gauges
    ! @param domain_time the time as a real number
    !
    SUBROUTINE write_current_time_series(point_gauges, domain_U, domain_time)
        CLASS(point_gauge_type), INTENT(INOUT):: point_gauges
        REAL(dp), INTENT(IN):: domain_U(:,:,:), domain_time
        INTEGER(ip) :: i
 
        CALL update_gauge_var(point_gauges, point_gauges%time_series_values, &
            domain_U, point_gauges%time_series_var)

#ifdef NONETCDF
        ! Write to home-brew binary format
        write(point_gauges%time_series_output_unit) &
            REAL(point_gauges%time_series_values, output_precision)
#else
        ! Write to ncdf
        ! We need to record how many time-steps are written
        point_gauges%netcdf_num_output_steps = point_gauges%netcdf_num_output_steps + 1
        
        ! Save the time
        call check(nf90_put_var(point_gauges%netcdf_gauge_output_file_ID, &
            point_gauges%netcdf_time_var_ID, domain_time, &
            start=[point_gauges%netcdf_num_output_steps]))

        ! Save all the time-series variables
        DO i = 1, size(point_gauges%time_series_var)
            call check(nf90_put_var(&
                point_gauges%netcdf_gauge_output_file_ID, &
                point_gauges%time_series_ncdf_iVar_ID(i), &
                point_gauges%time_series_values(:,i), &
                start=[1, point_gauges%netcdf_num_output_steps],&
                count=[point_gauges%n_gauges, 1]))
        END DO

        ! FIXME: Flushing the file is a crude way to ensure all data is written. If we can close the file
        ! when the simulation is finished, it shouldn't be needed.
        !if (flush_netcdf) call check(nf90_sync(point_gauges%netcdf_gauge_output_file_ID))
#endif

    END SUBROUTINE

    ! Initialise the gauges data structure and output file
    !
    ! Find the indices in domain_U associated with each gauge, given the
    ! lower left, dx, nx of the domain
    !
    ! @param point_gauges variable of type point_gauge_type
    ! @param domain_dx cell size [dx,dy] in the domain
    ! @param domain_nx number of cells [nx,ny] in the domain
    ! @param domain_U the array with dimensions [nx, ny, :] in which we look up values at the gauges
    ! @param netcdf_gauge_output_file name of output file
    ! @param attribute_names character vector of names for additional attributes to add to the netcdf file
    ! @param attribute_values character vector of values for the additional attributes in the netcdf file
    !
    SUBROUTINE initialise_gauges(point_gauges, domain_lower_left, domain_dx, &
        domain_nx, domain_U,&
        netcdf_gauge_output_file,&
        attribute_names, attribute_values)
        CHARACTER(charlen), INTENT(IN):: netcdf_gauge_output_file
        CLASS(point_gauge_type), INTENT(INOUT):: point_gauges
        REAL(dp), INTENT(IN):: domain_lower_left(:), domain_dx(:)
        INTEGER(ip), INTENT(IN):: domain_nx(:)
        REAL(dp), INTENT(IN):: domain_U(:,:,:)
        CHARACTER(charlen), OPTIONAL, INTENT(IN):: attribute_names(:), attribute_values(:)
    

        INTEGER(ip):: n_gauges, i, j, varind, xind, yind

        n_gauges = point_gauges%n_gauges
   
        ! Get indices of gauges on domain        
        DO i = 1, n_gauges
            point_gauges%site_index(:,i) = 1_ip + &
                floor((point_gauges%xy(:,i) - domain_lower_left)/domain_dx)

            ! Check it is inside the domain
            if((any(point_gauges%site_index(:,i) < 1_ip).OR. &
                any(point_gauges%site_index(:,i) > domain_nx))) then
                print*, 'Gauge coordinate ', point_gauges%xy(:,i), &
                    ' is outside the domain'
                call generic_stop()
            end if

        END DO
        
        ! Get static values -- could potentially write them here, and the
        ! deallocate the array, if conserving memory is very important
        CALL update_gauge_var(point_gauges, point_gauges%static_values, &
            domain_U, point_gauges%static_var)
        CALL update_gauge_var(point_gauges, point_gauges%time_series_values, &
            domain_U, point_gauges%time_series_var)

#ifndef NONETCDF
        ! Setup netcdf output
        if((present(attribute_names)).AND.(present(attribute_values))) then
            CALL setup_gauge_netcdf_output(point_gauges, netcdf_gauge_output_file,&
                attribute_names, attribute_values)
        else
            CALL setup_gauge_netcdf_output(point_gauges, netcdf_gauge_output_file)
        end if
#else
        ! Setup home-brew binary output
        !
        ! For simplicity, this routine only includes a netcdf filename as argument.
        ! If we are writing to a home-brew binary format, then we make file names from the 'netcdf' name
        point_gauges%static_output_file = TRIM(netcdf_gauge_output_file) // '_static' 
        point_gauges%time_series_output_file = TRIM(netcdf_gauge_output_file) // '_timeseries' 
        point_gauges%gauge_metadata_file = TRIM(netcdf_gauge_output_file) // '_metadata' 

        ! Save 'static' variables to home-brew binary
        open(newunit=point_gauges%static_output_unit, file=point_gauges%static_output_file, &
            access='stream', form='unformatted')
        DO i = 1, n_gauges
            write(point_gauges%static_output_unit) &
                REAL(point_gauges%xy(:,i), output_precision), &
                REAL(point_gauges%static_values(i,:), output_precision), &
                REAL(point_gauges%gauge_ids(i), output_precision) 
        END DO
        close(point_gauges%static_output_unit)

        ! Initialise time-varying output files
        open(newunit=point_gauges%time_series_output_unit, file=point_gauges%time_series_output_file, &
            access='stream', form='unformatted')

        ! Initialise gauge metadata file
        open(newunit=point_gauges%gauge_metadata_unit, file=point_gauges%gauge_metadata_file)
        write(point_gauges%gauge_metadata_unit, *) 'n_gauges : ', point_gauges%n_gauges
        write(point_gauges%gauge_metadata_unit, *) 'time_series_var: ', point_gauges%time_series_var
        write(point_gauges%gauge_metadata_unit, *) 'time_series_output_file: ', &
            TRIM(point_gauges%time_series_output_file)
        write(point_gauges%gauge_metadata_unit, *) 'static_var: ', point_gauges%static_var
        write(point_gauges%gauge_metadata_unit, *) 'static_output_file: ', &
            TRIM(point_gauges%static_output_file)
        close(point_gauges%gauge_metadata_unit)
#endif        
    END SUBROUTINE

    !
#ifndef NONETCDF
    !
    ! Convenience subroutine for netcdf error handling
    !
    SUBROUTINE check(error_status)
        INTEGER, INTENT ( IN) :: error_status
     
        if(error_status /= nf90_noerr) then
          print *, trim(nf90_strerror(error_status))
          call generic_stop()
        end if
    END SUBROUTINE

    !
    ! Create the gauge netcdf file and populate with key header information
    !
    ! @param point_gauges point_gauge_type
    ! @param netcdf_gauge_output_file filename where the gauges output will be stored
    !
    SUBROUTINE setup_gauge_netcdf_output(point_gauges, netcdf_gauge_output_file, &
        attribute_names, attribute_values)

        TYPE(point_gauge_type), INTENT(INOUT):: point_gauges 
        CHARACTER(charlen), INTENT(IN):: netcdf_gauge_output_file
        CHARACTER(charlen), OPTIONAL, INTENT(IN):: attribute_names(:)
        CHARACTER(charlen), OPTIONAL, INTENT(IN):: attribute_values(:)

        INTEGER:: iNcid, iDimStation_ID, iDimTime_ID, iDimLenStringName_ID, &
            iVarLON_ID, iVarLAT_ID, iVarTIME_ID, iVarGAUGEID_ID
        ! ID's for variables that we might store 'statically' (i.e. once at the
        ! start of the simulation)
        INTEGER:: iVar_STG_static_ID, iVar_UH_static_ID, iVar_VH_static_ID, &
            iVar_ELV_static_ID
        ! ID's for variables that we might store each output timestep 
        INTEGER:: iVar_STG_time_ID, iVar_UH_time_ID, iVar_VH_time_ID, &
            iVar_ELV_time_ID
        ! Up to 32 characters for site names (shorter than usual)
        ! Could set this based on input string name lengths 
        INTEGER:: iLenStringName = 32
        INTEGER(ip):: i, j, k
        CHARACTER(charlen):: err_mess

        point_gauges%netcdf_gauge_output_file = netcdf_gauge_output_file

        !! Using netcdf interface

        !! Create output file
        call check (nf90_create(netcdf_gauge_output_file, NF90_CLOBBER, iNcid))
    
        point_gauges%netcdf_gauge_output_file_ID = iNcid

        ! Define the dimensions. Try to follow CF standards for Orthogonal
        ! multidimensional array representation of time series
        call check(nf90_def_dim(iNcid, "station", point_gauges%n_gauges, iDimStation_ID))
        call check(nf90_def_dim(iNcid, "time", NF90_UNLIMITED, iDimTime_ID))
        call check(nf90_def_dim(iNcid, "NoCharStationName", iLenStringName, iDimLenStringName_ID))

        ! Define the variables

        ! Longitude
        call check( nf90_def_var(iNcid, "lon", NF90_REAL4, (/ iDimStation_ID /), iVarLON_ID) )
#ifdef SPHERICAL
        call check(nf90_put_att(iNcid, iVarLON_ID, "standard_name", &
            "longitude"))
        call check(nf90_put_att(iNcid, iVarLON_ID, "long_name", &
            "station_longitude"))
        call check(nf90_put_att(iNcid, iVarLON_ID, "units", "degrees_east"))
#else
        call check(nf90_put_att(iNcid, iVarLON_ID, "standard_name", &
            "projection_x_coordinate"))
        call check(nf90_put_att(iNcid, iVarLON_ID, "long_name", &
            "station_x_coordinate"))
        call check(nf90_put_att(iNcid, iVarLON_ID, "units", "m"))
#endif

        ! Latitude
        call check( nf90_def_var(iNcid, "lat", NF90_REAL4, (/ iDimStation_ID /), iVarLAT_ID) )
#ifdef SPHERICAL
        call check(nf90_put_att(iNcid, iVarLAT_ID, "standard_name", &
            "latitude"))
        call check(nf90_put_att(iNcid, iVarLAT_ID, "long_name", &
            "station_latitude"))
        call check(nf90_put_att(iNcid, iVarLAT_ID, "units", "degrees_north"))
#else
        call check(nf90_put_att(iNcid, iVarLAT_ID, "standard_name", &
            "projection_y_coordinate"))
        call check(nf90_put_att(iNcid, iVarLAT_ID, "long_name", &
            "station_y_coordinate"))
        call check(nf90_put_att(iNcid, iVarLAT_ID, "units", "m"))
#endif

        ! Time
        call check( nf90_def_var(iNcid, 'time', NF90_REAL4, (/ iDimTime_ID /), iVarTIME_ID) )
        point_gauges%netcdf_time_var_ID = iVarTIME_ID        
        call check(nf90_put_att(iNcid, iVarTIME_ID, "standard_name", &
            "time"))
        call check(nf90_put_att(iNcid, iVarTIME_ID, "long_name", &
            "time_from_start_of_simulation"))
        call check(nf90_put_att(iNcid, iVarTIME_ID, "units", "s"))

        ! GaugeID
        call check( nf90_def_var(iNcid, 'gaugeID', NF90_REAL4, (/iDimStation_ID/), iVarGAUGEID_ID))
        call check(nf90_put_att(iNcid, iVarGAUGEID_ID, "long_name", &
            "real_ID_for_each_station"))
        call check(nf90_put_att(iNcid, iVarGAUGEID_ID, "units", "-"))

        ! Define variables we store statically (i.e. once at the start)
        DO i = 1, size(point_gauges%static_var)
            j = point_gauges%static_var(i)

            ! Initial stage 
            if(j == STG) then
                call check( nf90_def_var(iNcid, "stage0", NF90_REAL4, [iDimStation_ID], iVar_STG_static_ID))
                call check(nf90_put_att(iNcid, iVar_STG_static_ID, "long_name", &
                    "water_surface_height_above_mean_sea_level_at_initial_condition"))
                call check(nf90_put_att(iNcid, iVar_STG_static_ID, "units", "m"))
            end if

            ! Initial UH
            if(j == UH) then
                call check( nf90_def_var(iNcid, "uh0", NF90_REAL4, [iDimStation_ID], iVar_UH_static_ID))
                call check(nf90_put_att(iNcid, iVar_UH_static_ID, "long_name", &
                    "x_velocity_multiplied_by_depth_at_initial_condition"))
                call check(nf90_put_att(iNcid, iVar_UH_static_ID, "units", "m^2/s"))
            end if

            ! Initial VH
            if(j == VH) then
                call check( nf90_def_var(iNcid, "vh0", NF90_REAL4, [iDimStation_ID], iVar_VH_static_ID))
                call check(nf90_put_att(iNcid, iVar_VH_static_ID, "long_name", &
                    "y_velocity_multiplied_by_depth_at_initial_condition"))
                call check(nf90_put_att(iNcid, iVar_VH_static_ID, "units", "m^2/s"))
            end if

            ! Initial elevation. Typically this is the only important one
            if(j == ELV) then
                call check( nf90_def_var(iNcid, "elevation0", NF90_REAL4, [iDimStation_ID], iVar_ELV_static_ID))
                call check(nf90_put_att(iNcid, iVar_ELV_static_ID, "long_name", &
                    "ground_level_altitude_above_mean_sea_level_at_initial_condition"))
                call check(nf90_put_att(iNcid, iVar_ELV_static_ID, "units", "m"))
            end if
        END DO

        ! Define variables we store each timestep.

        ! We need to store the netcdf Var ID's for these variables
        ALLOCATE(point_gauges%time_series_ncdf_iVar_ID(size(point_gauges%time_series_var)))

        DO i = 1, size(point_gauges%time_series_var)
            j = point_gauges%time_series_var(i)
            if (j == STG) then
                call check( nf90_def_var(iNcid, "stage", NF90_REAL4, [iDimStation_ID, iDimTime_ID], iVar_STG_time_ID))
                call check(nf90_put_att(iNcid, iVar_STG_time_ID, "long_name", &
                    "water_surface_height_above_mean_sea_level"))
                call check(nf90_put_att(iNcid, iVar_STG_time_ID, "units", "m"))
                point_gauges%time_series_ncdf_iVar_ID(i) = iVar_STG_time_ID
            end if
            if (j == UH) then
                call check( nf90_def_var(iNcid, "uh", NF90_REAL4, [iDimStation_ID, iDimTime_ID], iVar_UH_time_ID))
                call check(nf90_put_att(iNcid, iVar_UH_time_ID, "long_name", &
                    "x_velocity_multiplied_by_depth"))
                call check(nf90_put_att(iNcid, iVar_UH_time_ID, "units", "m^2/s"))
                point_gauges%time_series_ncdf_iVar_ID(i) = iVar_UH_time_ID
            end if
            if (j == VH) then
                call check( nf90_def_var(iNcid, "vh", NF90_REAL4, [iDimStation_ID, iDimTime_ID], iVar_VH_time_ID))
                call check(nf90_put_att(iNcid, iVar_VH_time_ID, "long_name", &
                    "y_velocity_multiplied_by_depth"))
                call check(nf90_put_att(iNcid, iVar_VH_time_ID, "units", "m^2/s"))
                point_gauges%time_series_ncdf_iVar_ID(i) = iVar_VH_time_ID
            end if
            if (j == ELV) then
                call check( nf90_def_var(iNcid, "elevation", NF90_REAL4, [iDimStation_ID, iDimTime_ID], iVar_ELV_time_ID))
                call check(nf90_put_att(iNcid, iVar_ELV_time_ID, "long_name", &
                    "ground_level_altitude_above_mean_sea_level"))
                call check(nf90_put_att(iNcid, iVar_ELV_time_ID, "units", "m"))
                point_gauges%time_series_ncdf_iVar_ID(i) = iVar_ELV_time_ID
            end if
        END DO

        ! This is the standard netcdf attribute for time-series
        call check(nf90_put_att(iNcid, nf90_global, "featureType", "timeSeries"))
   
        ! Here we add other attributes that might be useful (e.g. the name of the input stage and elevation) 
        if((present(attribute_names)).AND.(present(attribute_values))) then
            DO i = 1, size(attribute_names)
                call check(nf90_put_att(iNcid, nf90_global, attribute_names(i), attribute_values(i)))
            END DO
        end if

        ! Finish definitions so writing can begin
        call check(nf90_enddef(iNcid))
                
        ! Write lon/lat/gauge_ids
        call check(nf90_put_var(iNcid, iVarLON_ID , point_gauges%xy(1,:))) 
        call check(nf90_put_var(iNcid, iVarLAT_ID , point_gauges%xy(2,:))) 
        call check(nf90_put_var(iNcid, iVarGAUGEID_ID, point_gauges%gauge_ids(:)))

        ! Write static variables
        DO i = 1, size(point_gauges%static_var)
            j = point_gauges%static_var(i)
            if(j == STG) then
                call check(nf90_put_var(iNcid, iVar_STG_static_ID , point_gauges%static_values(:,i) )) 
            end if
            if(j == UH) then
                call check(nf90_put_var(iNcid, iVar_UH_static_ID , point_gauges%static_values(:,i) )) 
            end if
            if(j == VH) then
                call check(nf90_put_var(iNcid, iVar_VH_static_ID , point_gauges%static_values(:,i) ))
            end if
            if(j == ELV) then
                call check(nf90_put_var(iNcid, iVar_ELV_static_ID , point_gauges%static_values(:,i) )) 
            end if
        END DO

    END SUBROUTINE
#endif

    ! Cleanup point gauges
    ! We particularly use this to close netcdf files, which is essential to ensure
    ! they are fully written (unless we flush the file every timestep, which sounds inefficient)
    SUBROUTINE finalise_point_gauges(point_gauges)
        CLASS(point_gauge_type), INTENT(INOUT):: point_gauges

        if (allocated(point_gauges%xy)) then
#ifndef NONETCDF      
            ! Close the netcdf file 
            call check(nf90_close(point_gauges%netcdf_gauge_output_file_ID))
#else
            close(point_gauges%time_series_output_unit)
#endif
        end if

    END SUBROUTINE

    !
    ! Convenienve routine for testing
    !
    SUBROUTINE assert_test(test)
        LOGICAL, INTENT(IN):: test
        if(test) then
            print*, 'PASS'
        else
            print*, 'FAIL'
        end if
    END SUBROUTINE


    SUBROUTINE test_point_gauge_mod
        TYPE(point_gauge_type):: point_gauges
        REAL(dp):: domain_ll(2), domain_dx(2)
        INTEGER(ip):: domain_nx(2), time_series_var_indices(3), static_var_indices(1)
        INTEGER(dp):: domain_nvar
        REAL(dp):: xy_coords(2,6)
        REAL(dp), ALLOCATABLE:: domain_U(:,:,:)
        REAL(dp):: gauge_ids(6)

        INTEGER(ip) :: i, j, k
        CHARACTER(charlen), PARAMETER:: netcdf_gauges_file = 'test_gauge_output.nc'
        CHARACTER(charlen):: attribute_names(2), attribute_values(2)

        attribute_names(1) = 'unit_test_attr1'
        attribute_names(2) = 'unit_test_attr2'
        attribute_values(1) = 'unit_test_values1'
        attribute_values(2) = 'unit_test_values2'

        ! Make a skeleton of a domain
        domain_nvar = 4_ip
        domain_nx = [104_ip, 104_ip]
        domain_dx = [1.0_dp, 2.0_dp]
        domain_ll = [-100.0_dp, 250.0_dp]
        ALLOCATE(domain_U(domain_nx(1), domain_nx(2), domain_nvar))

        ! Make some gauge coordinates
        xy_coords(:,1) = [0._dp, 0._dp] + domain_ll
        xy_coords(:,2) = [50.0_dp, 20.0_dp] + domain_ll
        xy_coords(:,3) = [20.0_dp, 50.0_dp] + domain_ll
        xy_coords(:,4) = [20.0_dp, 200.0_dp] + domain_ll
        xy_coords(:,5) = [90.0_dp, 10.0_dp] + domain_ll
        xy_coords(:,6) = [90.0_dp, 10.0_dp] + domain_ll

        gauge_ids = (/ (i, i=1,6) /)

        ! Make up some domain variables that are easy to id
        DO k = 1, domain_nvar
            DO j = 1, domain_nx(2)
                DO i = 1, domain_nx(1)
                    domain_U(i, j, k) = i + j**2 + k**3
                END DO
            END DO
        END DO
      
        ! Choose some indices to store statically/dynamically 
        time_series_var_indices = [1_ip, 2_ip, 4_ip]
        static_var_indices = 3_ip 

        call point_gauges%allocate_gauges(xy_coords, time_series_var_indices, &
            static_var_indices, gauge_ids)

        !! Tests
        call assert_test(all(point_gauges%gauge_ids == gauge_ids))
        call assert_test(all(point_gauges%xy == xy_coords))
        call assert_test(point_gauges%n_gauges == size(xy_coords(1,:)))
        call assert_test(size(point_gauges%static_values(:,1)) == 6)
        call assert_test(size(point_gauges%static_values(1,:)) == 1)
        call assert_test(size(point_gauges%time_series_values(1,:)) == 3)
        call assert_test(size(point_gauges%time_series_values(:,1)) == 6)
        call assert_test(all(point_gauges%static_var == static_var_indices))
        call assert_test(all(point_gauges%time_series_var == time_series_var_indices))
        call assert_test(size(point_gauges%site_index(:,1)) == 2)
        call assert_test(size(point_gauges%site_index(1,:)) == 6)

        ! Check that gauges can be initialised
        call point_gauges%initialise_gauges(domain_ll, domain_dx, domain_nx, domain_U, &
            netcdf_gauges_file, attribute_names=attribute_names, attribute_values=attribute_values)
        call point_gauges%write_current_time_series(domain_U, 0.0_dp)

        call assert_test(all(point_gauges%site_index(1,:) == &
            floor((point_gauges%xy(1,:) - domain_ll(1))/domain_dx(1)) + 1))
        call assert_test(all(point_gauges%site_index(2,:) == &
            floor((point_gauges%xy(2,:) - domain_ll(2))/domain_dx(2)) + 1))

        call assert_test(all(point_gauges%site_index <= 104))
        call assert_test(all(point_gauges%site_index > 0))

        call assert_test(all(point_gauges%static_values(:,1) == &
            point_gauges%site_index(1,:) + point_gauges%site_index(2,:)**2 + 3**3))
       
        call assert_test(all(point_gauges%time_series_values(:,1) == &
            point_gauges%site_index(1,:) + point_gauges%site_index(2,:)**2 + 1**3))
        call assert_test(all(point_gauges%time_series_values(:,2) == &
            point_gauges%site_index(1,:) + point_gauges%site_index(2,:)**2 + 2**3))
        call assert_test(all(point_gauges%time_series_values(:,3) == &
            point_gauges%site_index(1,:) + point_gauges%site_index(2,:)**2 + 4**3))

        ! Check the update
        domain_U = domain_U - 5.0_dp
        call point_gauges%update_gauge_var(point_gauges%time_series_values, &
            domain_U, point_gauges%time_series_var)
        call point_gauges%write_current_time_series(domain_U, 1.5_dp)

        call assert_test(all(point_gauges%time_series_values(:,1) == &
            point_gauges%site_index(1,:) + point_gauges%site_index(2,:)**2 + 1**3 - 5.0_dp))
        call assert_test(all(point_gauges%time_series_values(:,2) == &
            point_gauges%site_index(1,:) + point_gauges%site_index(2,:)**2 + 2**3 - 5.0_dp))
        call assert_test(all(point_gauges%time_series_values(:,3) == &
            point_gauges%site_index(1,:) + point_gauges%site_index(2,:)**2 + 4**3 - 5.0_dp))

        call point_gauges%finalise()

    END SUBROUTINE

END MODULE
