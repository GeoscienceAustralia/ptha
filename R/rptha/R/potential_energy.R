#' Available potential energy of the sea surface.
#'
#' The available potential energy is computed as the following integral in space: \cr
#'     integral [ 0.5 * gravity * seawater_density * (sea_surface - MSL)^2 ] dA \cr
#' where dA is the 2D area element. \cr
#' This function takes a matrix giving the sea_surface, as well as the cell areas (which should
#' be a matrix unless the cell areas are constant). You must ensure that the sea_surface is NA
#' in dry areas.
#'
#' @param sea_surface matrix giving the initial sea state. IT SHOULD HAVE NA VALUES IN DRY AREAS.
#' @param cell_areas_m2 The cell area associated with each entry of sea_surface
#' in meters**2. This can be a constant (e.g. for regular cells in Cartesian coordinates) or
#' it can be a matrix of the same size as sea_surface
#' @param gravity Gravitational acceleration at the earth's surface (m/s^2)
#' @param seawater_density Density of seawater in kg/m^3
#' @param MSL the mean sea level
#' @return The available potential energy in units of joules (kg m^2/s^2)
#' @export
#' @examples
#' # Read a raster initial condition -- this is in lonlat coordinates, and is 
#' # already NA on land. It was created from this paper: https://doi.org/10.1029/2018JB016996
#' raster_file = system.file('extdata/Ho_Chile1960_initial_displacement_NA_on_land.tif', 
#'     package='rptha')
#' chile1960_free_surface_NA_on_land = raster(raster_file)
#' # Get the cell areas using the raster package's "area" function. This correctly
#' # recognizes that the raster is lonlat
#' cell_areas_m2 = area(chile1960_free_surface_NA_on_land) * 1e+03 * 1e+03 # convert km^2 to m^2
#' # Compute the available potential energy. It should be about 6.8x10^15
#' ape = sea_surface_available_potential_energy(
#'     as.matrix(chile1960_free_surface_NA_on_land),
#'     as.matrix(cell_areas_m2),
#'     gravity=9.8, seawater_density=1024, MSL=0)
#' stopifnot(signif(ape, 2) == 6.8e+15)
#' 
sea_surface_available_potential_energy<-function(
    sea_surface, 
    cell_areas_m2, 
    gravity=9.8, 
    seawater_density=1024, 
    MSL=0){

    if(!(length(dim(sea_surface))==2)){
        stop('dim(sea_surface) should have a length of 2')
    }

    if(!(length(dim(cell_areas_m2))==0)){
        stopifnot(all(dim(cell_areas_m2) == dim(sea_surface)))
    }

    ape = sum(gravity/2 * seawater_density * cell_areas_m2 * 
              (sea_surface-MSL)**2, na.rm=TRUE)

    return(ape)

}
