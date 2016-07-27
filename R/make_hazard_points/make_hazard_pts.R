#
# This script can be used to make hazard points for PTHA from a DEM.
#
# It requires user editing of the input_data below, and potentially also
# creation of some of the input files (e.g. no_clip_zone or haz_pts_mask)
#
# It also requires that a range of standard gdal tools are in your path and
# can be called linux-style.
#
# The key output is a file of hazard point locations
#
# To run this script interactively and check progress, start R and run:
# source('make_hazard_pts.R', echo=TRUE, max.deparse.length=Inf)
#
# For non-interactive usage, use the following command (from the commandline)
# Rscript make_hazard_pts.R


###############################################################################
## INPUT DATA 

# Main raster file used for contour creation
raster_infile = '../../../../DATA/ELEV/GEBCO_08/gebco_08_W-180-E180-S-90-N90-1mx1m.nc'

# Hazard contour level (m). Elevation at which we 'desire' hazard points to be
# located (e.g. -100 for 100m water depth). 
hazard_contour_level = -100 

# Spacing of hazard points along the contour (km)
hazard_pt_spacing = 25 

# Elevation defining a contour called 'the coast'
# This value will place an upper bound on the hazard point elevation (although
# usually the latter will be closer to hazard_contour_level)
coast_contour_level = -0.001 

# We will reject hazard points that are within this many raster cells of 'the
# coast'
coast_buffer_ncell = 3 

# We will accept hazard points that are within this many raster cells from 'the
# coast'. 
hazard_buffer_ncell = 44 

# Polygon shapefile name defining the region where we don't clip small islands
# Set to NULL if you don't want to use this
no_clip_zone = '../../../../DATA/ELEV/ISLAND_CLIP_LAYER/ISLAND_CLIP_LAYER.shp' 

# Outside of the no-clip-zone, a coast contour must enclose at least
# coast_contour_removal_area_threshold (km^2), otherwise we remove it on
# creation (so it will not be surrounded by hazard points). This is to avoid
# having too many hazard points associated with tiny islands.
coast_contour_removal_area_threshold = 100 

# 'Land' values are set to this value before hazard contour computation. It
# helps the contour algorithm avoid crossing land values, and probably does not
# need to be changed
land_value = 10000 

# Name of a line shapefile where we will have extra hazard points sampled.
# We use this in locations where our algorithm does not produce hazard points,
# even though we want them. Set to NULL to not use anything.
extra_manual_haz_lines = NULL

# Name of a polygon shapefile, used as a mask. We will not include hazard
# points inside this polygon. Set to NULL if there is no mask
haz_pts_mask = '../../../../DATA/ELEV/HAZ_PT_REMOVAL_REGION/HAZ_PT_REMOVAL_REGION.shp' 

# FINAL STEP: Translate hazard points so this is the lower left longitude -- to
# match with the tsunami model domain
haz_pt_lowerleft = -45 

## END INPUT DATA
###############################################################################
#
# Get libraries
library(rptha)
# Utility routines for dealing with SpatialLines
cu = new.env()
sys.source('contour_util.R', cu)
# Utility routines for dealing with SpatialPoints
pu = new.env()
sys.source('point_util.R', pu)

##############################################################################
#
# Read raster, reformat to Geotif for GDAL operations, and make a zero contour

dem = raster(raster_infile)

# Use gdal to convert to geotif -- assume WGS84
#out_tif_dir = paste0('./OUTPUTS/', basename(dirname(raster_infile)))
out_tif_dir = './OUTPUTS/raster_tmp'
dir.create(out_tif_dir, showWarnings=FALSE, recursive=TRUE)
temp_tif = paste0(out_tif_dir, '/', 
    strsplit(basename(raster_infile), split='\\.')[[1]][1], '.tif')
gdal_translate_command = paste('gdal_translate -a_srs EPSG:4326', raster_infile, 
    temp_tif, ' -co COMPRESS=DEFLATE')
system(gdal_translate_command)

# Read in the (properly georeferenced) data, get zero contour
dem = raster(temp_tif)
dem_0m = cu$gdal_contour(temp_tif, contour_levels=-0.0001, 
    out_dsn='OUTPUTS/ZERO_CONTOUR', contour_shp_name='ZERO_CONTOUR')
###############################################################################


###############################################################################
#
# Cut contourLines with length < coast_contour_removal_length_threshold.
# Except in the area inside the 'no-clip-zone', which protects e.g. small
# pacific countries
#
dem_0m_poly = cu$SpatialLinesDF2Polygons(dem_0m) 

# NOTE: Be careful of self-intersections for areaPolygon -- can give crazy
# answers in rare situations (e.g. antarctica)
dem_0m_polyArea = geosphere::areaPolygon(dem_0m_poly)/1e+06 

# Find a point inside each contour polygon -- so we can check if it is in the
# no-clip-zone
if(!is.null(no_clip_zone)){
    no_clip_area = readOGR(dsn=no_clip_zone, 
        layer=gsub('.shp', '', basename(no_clip_zone)))
    no_clipping_pt = SpatialPoints(coordinates(dem_0m_poly), 
        proj4string=CRS(proj4string(no_clip_area)))
    in_noclip_zone = !is.na(over(no_clipping_pt, no_clip_area)) 
}else{
    # There are no points 'protected' from removal
    in_noclip_zone = dem_0m_polyArea * 0
}

# Contours to keep
keep_islands = which((dem_0m_polyArea > coast_contour_removal_area_threshold ) |
    (in_noclip_zone == 1))

dem_0m_trim = dem_0m[keep_islands,]

writeOGR(dem_0m_trim, dsn='OUTPUTS/ZERO_CONTOUR_TRIMMED', 
    layer='ZERO_CONTOUR_TRIMMED', driver='ESRI Shapefile',
    overwrite=TRUE)

# Sample points along the coastal contour -- we will later use these to
# identify hazard_lines that do not loop around land
dem_0m_trim_pts = rptha::approxSpatialLines(dem_0m_trim, 
    spacing=hazard_pt_spacing,longlat=TRUE)

###############################################################################

###############################################################################
#
# Create a raster that has 1's within a 'band' where we would like all out
# hazard points to be (i.e. within a certain distance from the coast), and
# zeros elsewhere. This is a many step process
#

# Burn the zero contour into a raster. Initialise its values to 0.0, and
# define these values to be NA (so that gdal_fillnodata.py can use them).
gdal_rasterize_command = paste(
    "gdal_rasterize -burn 1.0 -at -init 0.0",
    " -tr", res(dem)[1], res(dem)[2], 
    '-te ', extent(dem)@xmin, extent(dem)@ymin,
            extent(dem)@xmax, extent(dem)@ymax,
    '-a_nodata 0.0',
    ' OUTPUTS/ZERO_CONTOUR_TRIMMED/ZERO_CONTOUR_TRIMMED.shp',
    'OUTPUTS/raster_tmp/raster_tmp.tif',
    ' -co COMPRESS=DEFLATE -co BIGTIFF=YES')
system(gdal_rasterize_command)

# Make a copy of this raster -- we will edit 2 x independently
system('cp OUTPUTS/raster_tmp/raster_tmp.tif OUTPUTS/raster_tmp/raster_tmp3.tif')

# Use gdal_fillnodata to expand the zero contour into a region
gdal_fillnodata_command = paste(
    'gdal_fillnodata.py -md ', hazard_buffer_ncell, 
    'OUTPUTS/raster_tmp/raster_tmp.tif',
    ' -co COMPRESS=DEFALTE -co BIGTIFF=YES')
system(gdal_fillnodata_command)

# Change the nodata value in the tif to negative 1.0 -- so that we can use
# gdal_calc with the 0 and 1 values treated as numbers
gdal_changenodata_command = paste(
    'gdal_translate -a_nodata -1.0',
    ' OUTPUTS/raster_tmp/raster_tmp.tif',
    ' OUTPUTS/raster_tmp/raster_tmp2.tif',
    ' -co COMPRESS=DEFLATE -co BIGTIFF=YES')
system(gdal_changenodata_command)

# Cut-out points that are TOO NEAR the zero contour. This will be entirely
# contained within the previously created hazard points. It is best not to have
# hazard points very close to the coast.
gdal_fillnodata_command = paste(
    'gdal_fillnodata.py -md ', coast_buffer_ncell, 
    'OUTPUTS/raster_tmp/raster_tmp3.tif',
    ' -co COMPRESS=DEFALTE -co BIGTIFF=YES')
system(gdal_fillnodata_command)

# Fix nodata values
gdal_changenodata_command = paste(
    'gdal_translate -a_nodata -1.0',
    ' OUTPUTS/raster_tmp/raster_tmp3.tif',
    ' OUTPUTS/raster_tmp/raster_tmp4.tif',
    ' -co COMPRESS=DEFLATE -co BIGTIFF=YES')
system(gdal_changenodata_command)

gdal_calc_command=paste(
    'gdal_calc.py -A OUTPUTS/raster_tmp/raster_tmp2.tif',
    ' -B OUTPUTS/raster_tmp/raster_tmp4.tif',
    " --calc='(A-B)'",
    ' --outfile=OUTPUTS/raster_tmp/raster_tmp5.tif',
    ' --co COMPRESS=DEFLATE BIGTIFF=YES')
system(gdal_calc_command)

###############################################################################

###############################################################################
#
# Make a raster, which we can contour to get our hazard points in the right
# place
#
# A = Elevation data;
# B = raster with 1.0 in the 'areas we want hazard points', and 0.0 elsewhere
# C = raster with 1.0 in the 'near coast areas where we don't want hazard
#     points', and 0 elsewhere
#
# What we want: If (B==1) -- we want A
#               If (B==0) and A<0 -- We want -land_thresh*(a_bit_less_than_1)
#               If (B==0) and A>0 OR C==1 -- Wa want land_thresh (or another very large number)

land_thresh = land_value
outfile_gdc = paste0(dirname(temp_tif), '/gdc_proc.tif')
gdal_calc_command = paste(
    "gdal_calc.py -A", temp_tif, 
    " -B OUTPUTS/raster_tmp/raster_tmp5.tif", 
    " -C OUTPUTS/raster_tmp/raster_tmp4.tif",
    " --calc=' A*(B==1) + ", 
    land_thresh, "*(((C==1) + (A>", coast_contour_level,"))>0) + ",
    -land_thresh*0.9, "*((B+C)==0)*(A<=", coast_contour_level, ")'",
    " --outfile", outfile_gdc, 
    " --co=COMPRESS=DEFLATE")
system(gdal_calc_command)

###############################################################################

###############################################################################
# Make hazard points
# Make a contour of the new DEM
contour_depth = hazard_contour_level
dem_haz_cont = cu$gdal_contour(outfile_gdc, contour_levels=contour_depth, 
    out_dsn='OUTPUTS/HAZARD_CONTOUR', 
    contour_shp_name='HAZARD_CONTOUR')

# Convert dem_has_cont to polygons, and remove areas not containing coastal
# points
dem_haz_poly = cu$SpatialLinesDF2Polygons(dem_haz_cont)
# Find those with coastal points inside
FirstPointIndex = over(dem_haz_poly, dem_0m_trim_pts)
dem_haz_poly_clip = dem_haz_poly[!is.na(FirstPointIndex[[1]]),]
dem_haz_lines = as(dem_haz_poly_clip,'SpatialLines')

# Combine with manual lines to deal with areas that our algorithm 'undesirably'
# removes.
if(!is.null(extra_manual_haz_lines)){
    extra_haz_lines = readOGR(extra_manual_haz_lines, 
        layer=gsub('.shp', '', basename(extra_manual_haz_lines)))
    extra_haz_lines = as(extra_haz_lines, 'SpatialLines')
    dem_haz_lines = rbind(dem_haz_lines, extra_haz_lines)
}

## Get points spaced at 'hazard_pt_spacing' km
haz_pts = rptha::approxSpatialLines(dem_haz_lines, spacing=hazard_pt_spacing, 
    longlat=TRUE)

haz_pts@data = cbind(haz_pts@data, data.frame('elev'=dem[haz_pts]))

# Ensure all points have elev < 0 and non NA elev.
# At this stage, we will still have occasional points with terrestrial
# elevation. This can be due to e.g. sites where we jump from 'deep ocean' to
# 'land' [which becomes more likely due to our island removal]. The contouring
# algorithm can interact with this to produce a few points with positive
# elevation. E.G. in one instance I had ~ 18000 hazard points, with 79 being
# terrestrial.
# Also, NA hazard points can happen at the raster boundary, remove them
haz_pts = haz_pts[haz_pts$elev < 0 & is.finite(haz_pts$elev),]

# Write out
writeOGR(haz_pts, dsn='OUTPUTS/HAZ_PTS', layer='HAZ_PTS',
         driver='ESRI Shapefile', overwrite=TRUE)

# Edit the hazard points to A) remove from mask polygon, and B) adjust
# lower-left coordinate
haz_new = pu$cut_hazpts_in_poly(haz_pts, haz_pts_mask, 
    'HAZ_NEW', lower_left=haz_pt_lowerleft)

pu$haz_pts_2_ursga_format(haz_new, 
    outfile=paste0('OUTPUTS/haz_pts_W', haz_pt_lowerleft,'.txt'))

