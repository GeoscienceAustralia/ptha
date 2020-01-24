#
# This script shows how to extract 3D displacement vectors at a given point
# (target_pt) for a set of ptha earthquake events
#

# We want to get 3D displacements at this point (lon/lat)
target_pt = c(360-175.1982, -21.1790)
source_zone = 'kermadectonga2'

# Get scripts to read PTHA results, just inside the 'ptha' repository under
# ptha_access
source('/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R', chdir=TRUE)

# Get Kermadec-Tonga2 earthquake scenarios in PTHA.
kt2 = get_source_zone_events_data(source_zone=source_zone)


# The unit-source data, stored in tif format
vertical_disp_files = Sys.glob(paste0('OUTPUTS/Unit_source_data/', source_zone, '/', source_zone, '_vertical_displacement*.tif'))
northing_disp_files = Sys.glob(paste0('OUTPUTS/Unit_source_data/', source_zone, '/', source_zone, '_northing_displacement*.tif'))
easting_disp_files = Sys.glob(paste0('OUTPUTS/Unit_source_data/', source_zone, '/', source_zone, '_easting_displacement*.tif'))

# Function to get the unit-source results at our target location
get_displacements_due_to_unit_source<-function(
    target_pt,
    unit_source_alongstrike_index,
    unit_source_downdip_index){

    # Find the files for this unit source
    vertical_rast = grep(
        paste0('_vertical_displacement_', unit_source_downdip_index, '_', unit_source_alongstrike_index, '.tif'),
        vertical_disp_files)
    stopifnot(length(vertical_rast) == 1)
    vertical_rast = vertical_disp_files[vertical_rast]

    easting_rast = grep(
        paste0('_easting_displacement_', unit_source_downdip_index, '_', unit_source_alongstrike_index, '.tif'),
        easting_disp_files)
    stopifnot(length(easting_rast) == 1)
    easting_rast = easting_disp_files[easting_rast]

    northing_rast = grep(
        paste0('_northing_displacement_', unit_source_downdip_index, '_', unit_source_alongstrike_index, '.tif'),
        northing_disp_files)
    stopifnot(length(northing_rast) == 1)
    northing_rast = northing_disp_files[northing_rast]

    # Get the x/y/z values
    xyz=rep(NA,3)
    xyz[1] = extract(raster(easting_rast), matrix(target_pt, nrow=1))
    xyz[2] = extract(raster(northing_rast), matrix(target_pt, nrow=1))
    xyz[3] = extract(raster(vertical_rast), matrix(target_pt, nrow=1))

    return(xyz)

}    

# Make a Nx3 matrix containing x/y/z displacements for the N unit_sources
xyz_displacement_unit_sources = matrix(NA, ncol=3, nrow=nrow(kt2$unit_source_statistics))
for(i in 1:nrow(xyz_displacement_unit_sources)){
    #print(i)
    xyz_displacement_unit_sources[i,1:3] = get_displacements_due_to_unit_source(
        target_pt,
        kt2$unit_source_statistics$alongstrike_number[i],
        kt2$unit_source_statistics$downdip_number[i])
    # Many unit-sources produce a zero displacement, which reflects that in PTHA18
    # we limit the region of the displacement calculations (depending on the
    # unit-source depth as well as horizontal distances). This is reasonable because
    # Okada displacements are concentrated near the source, albeit more spread out with deeper earthquakes.
    # The truncation options can be changed, see config.R.
}

# Now we can compute the displacement for all earthquake event
xyz_displacement_events = matrix(NA, ncol=3, nrow=nrow(kt2$events))
for(i in 1:nrow(kt2$events)){

    # Extract the slip from the formatted string (has slip separated by '_')
    event_slip = as.numeric(strsplit(kt2$events$event_slip_string[i], '_')[[1]])

    # Extract the unit-source indices from the formatted string (integers separated by '-')
    # Each integer corresponds to a row-index in kt2$unit_source_statistics, and thus a
    # row of xyz_displacement_unit_sources
    event_unit_source_inds = as.numeric(strsplit(kt2$events$event_index_string[i], '-')[[1]])

    stopifnot(length(event_slip) == length(event_unit_source_inds))

    # Linear summation to get the displacements
    xyz_displacement_events[i,1] = sum(xyz_displacement_unit_sources[event_unit_source_inds,1]*event_slip)
    xyz_displacement_events[i,2] = sum(xyz_displacement_unit_sources[event_unit_source_inds,2]*event_slip)
    xyz_displacement_events[i,3] = sum(xyz_displacement_unit_sources[event_unit_source_inds,3]*event_slip)
}

save.image('3D_displacements_R_image.Rdata')
