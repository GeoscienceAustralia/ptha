#
# Make a single netcdf file, containing the 'peak stage' at all gauges due to
# each individual unit source. This makes it easy/fast to access the data e.g.
# for plotting.
#
library(rptha)

output_dir = '/g/data/fj6/PTHA/AustPTHA_1/EVENT_RATES'

# All unit source files
all_uss_files = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/unit_source_statistics*.nc')
all_uss = sapply(all_uss_files, read_table_from_netcdf, simplify=FALSE)
# Put in a single data.frame
uss = do.call(rbind, all_uss)
# Add a sourcename column
uss$sourcename = basename(dirname(dirname(dirname(dirname(dirname(uss$tide_gauge_file))))))

# Read lon/lat/elev/gaugeID
get_lon_lat_elev_id<-function(tide_gauge_file){
    fid = nc_open(tide_gauge_file)
    lon = ncvar_get(fid, 'lon')
    lat = ncvar_get(fid, 'lat')
    elev = ncvar_get(fid, 'elevation0')
    gaugeID = ncvar_get(fid, 'gaugeID')
    nc_close(fid)

    output = data.frame(lon=as.numeric(lon), lat=as.numeric(lat), 
        elev=as.numeric(elev), gaugeID=as.numeric(gaugeID))


    return(output)
}

# Get coordinates for all sites
site_info = get_lon_lat_elev_id(uss$tide_gauge_file[1])

# Function to get peak stages at all hazard points, from a single unit-source
get_unit_source_peak_stage<-function(i){
    tgf = uss$tide_gauge_file[i]

    fid = nc_open(tgf)
    stage_block = ncvar_get(fid, 'stage')
    nc_close(fid)

    stage_max = apply(stage_block, 2, max)    
    rm(stage_block); gc()
    #peak_stage[,i] = stage_max
    return(stage_max)
}

#
# Read the peak stages
#
library(parallel)
cl = makeForkCluster(nnodes=16)
peak_stages = parLapply(cl, as.list(1:length(uss[,1])), get_unit_source_peak_stage)
new_peak_stages = do.call(cbind, peak_stages)
stopCluster(cl)

#
# Save to netcdf
#

# Get a dimension for the station
fid_tg = nc_open(uss$tide_gauge_file[1])
station_dim = fid_tg$dim$station
nc_close(fid_tg)

# Get a dimension for the unit-sources
unit_source_dim = ncdim_def('unitsource_index', "", 1:length(uss[,1]), unlim=FALSE)
# Get a 'maximum characters' dim
char_limit = max(max(nchar(uss$tide_gauge_file)), max(nchar(uss$initial_condition_file)))
max_nchar_dim = ncdim_def('max_nchar', '', 1:char_limit, unlim=FALSE)

all_nc_vars = list()

#
# Make numeric unit-source variables
#
var_names = c('lon_c', 'lat_c', 'depth', 'strike', 'dip', 'rake', 'slip', 
    'length', 'width', 'downdip_number', 'alongstrike_number', 
    'subfault_number', 'max_depth')
var_units = c('degreesE', 'degreesN','km', 'degrees', 'degrees', 'degrees', 
    'm', 'km', 'km', '', '', '', 'km')
var_longnames = c('unit source longitude', 
    'unit source latitude', 
    'unit source centroid depth (approximate only)', 
    'unit source strike (approximate only)', 
    'unit source dip (approximate only)', 
    'unit source rake (always -90 or 90)',
    'unit source slip (always 1 m)',
    'unit source length (approximate only)',
    'unit source width down-dip (defined so that length x width = area)',
    'a downdip index on the source-zone',
    'an along-strike index on the source-zone',
    'a unit-source index',
    'maximum depth on the unit-source')
for(i in 1:length(var_names)){
    all_nc_vars[[var_names[i]]] = ncvar_def(
        var_names[i],
        var_units[i],
        dim=unit_source_dim,
        longname=var_longnames[i])
}

#
# Append name variables
#
var_names = c('sourcename', 'initial_condition_file', 'tide_gauge_file')
var_units = c('', '', '')
var_longnames = c('', '', '')
for(i in 1:length(var_names)){
    # 
    all_nc_vars[[ var_names[i] ]] = ncvar_def(
        var_names[i],
        var_units[i],
        dim=list(max_nchar_dim, unit_source_dim),
        longname=var_longnames[i],
        prec='char')
}

#
# Append tide-gauge-variables
#
var_names = c('lon', 'lat', 'elev', 'gaugeID')
var_units = c('degreesE', 'degreesN', 'm', '')
var_longnames = c('gauge longitude', 'gauge latitude', 'gauge elevation above MSL', 'gauge ID')
for(i in 1:length(var_names)){

    all_nc_vars[[var_names[i]]] = ncvar_def(
            var_names[i],
            var_units[i],
            dim=list(station_dim),
            longname=var_longnames[i],
            prec='float')

}


# 
# Append the peak stage variable
#
all_nc_vars[['max_stage']] = ncvar_def('max_stage', 'm', dim=list(unit_source_dim, station_dim), 
        longname='maximum stage from the unit-source', prec='float')

# Make the file
fid = nc_create(paste0(output_dir, '/all_unit_source_wave_heights.nc'), all_nc_vars)

# Add the unit-source variables
for(i in 1:length(names(uss))){
    ncvar_put(fid, names(uss)[i], uss[,i])
}
# Add the station variables
for(i in 1:length(names(site_info))){
    ncvar_put(fid, all_nc_vars[[ names(site_info)[i] ]], site_info[,i], start=1, count=length(site_info[,1]))
}
# Add the max stage
ncvar_put( fid, 'max_stage', t(new_peak_stages), start=c(1,1), count=c(ncol(new_peak_stages), nrow(new_peak_stages)) )

nc_close(fid)


