#
# Single station stage exceedance rate computation
#
# This is actually a slow computational method ('quick' means 'quick to code'), but is useful
# to cross check the other results

lon = 151.41
lat = -34.08

tsunami_files = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_tsunami_*.nc')

library(rptha)

# Get hazard points -- faster to not use the '_tsunami' file
fid = nc_open(gsub('_tsunami', '', tsunami_files[1], fixed=TRUE), readunlim=FALSE)
n = length(fid$dim$station$vals)
hp = data.frame(
    lon     = rep(NA, n),
    lat     = rep(NA, n),
    elev    = rep(NA, n),
    gaugeID = rep(NA, n)
    )
for(i in 1:n){
    hp$lon[i] = ncvar_get(fid, 'lon', start=i, count=1)
    hp$lat[i] = ncvar_get(fid, 'lat', start=i, count=1)
    hp$elev[i] = ncvar_get(fid, 'elev', start=i, count=1)
    hp$gaugeID[i] = ncvar_get(fid, 'gaugeID', start=i, count=1)
    if(i%%100 == 1) print(i)
}
nc_close(fid)

# Find index of point nearest to lon/lat
ni = lonlat_nearest_neighbours(cbind(lon, lat), cbind(hp$lon, hp$lat))

# Get stage and rates, for each source

stage_rate = list()
names(stage_rate) = basename(dirname(dirname(tsunami_files)))
for(i in 1:length(tsunami_files)){
    print(i)
    fid = nc_open(tsunami_files[i], readunlim=FALSE)

    event_rate = ncvar_get(fid, 'event_rate_annual')
    event_rate_upper = ncvar_get(fid, 'event_rate_annual_upper_ci')
    event_rate_lower = ncvar_get(fid, 'event_rate_annual_lower_ci')
    peak_stage = ncvar_get(fid, 'max_stage', start=c(1, ni), count=c(-1,1))
    site = rep(basename(dirname(dirname(tsunami_files[i]))), length=length(event_rate))

    stage_rate[[i]] = data.frame(
        event_rate = event_rate,
        event_rate_upper = event_rate_upper,
        event_rate_lower = event_rate_lower,
        peak_stage = peak_stage,
        site = site)

    nc_close(fid)
}

# Back-calculate the stage-vs-rate curves
stage_rate_all = do.call(rbind, stage_rate)

odr = rev(order(stage_rate_all$peak_stage))

stg = stage_rate_all$peak_stage[odr]
er = cumsum(stage_rate_all$event_rate[odr])
er_up = cumsum(stage_rate_all$event_rate_upper[odr])
er_lo = cumsum(stage_rate_all$event_rate_lower[odr])

#
# Plot the data
#
plot(stg, er, t='l', log='xy', xlim=c(0.01, max(stg)))
grid()
points(stg, er_up, t='l', col='red')
points(stg, er_lo, t='l', col='red')

# Compare with the values in the file
fid = nc_open('tsunami_stage_exceedance_rates_sum_over_all_source_zones.nc', 
    readunlim=FALSE)
stages = fid$dim$stage$vals
ers = ncvar_get(fid, 'stochastic_slip_rate'            , start=c(1,ni), count=c(-1,1))
ers_up = ncvar_get(fid, 'stochastic_slip_rate_upper_ci', start=c(1,ni), count=c(-1,1))
ers_lo = ncvar_get(fid, 'stochastic_slip_rate_lower_ci', start=c(1,ni), count=c(-1,1))
nc_close(fid)

points(stages, ers, pch=19, cex=0.6, col='brown')
points(stages, ers_up, pch=19, cex=0.6, col='pink')
points(stages, ers_lo, pch=19, cex=0.6, col='pink')


