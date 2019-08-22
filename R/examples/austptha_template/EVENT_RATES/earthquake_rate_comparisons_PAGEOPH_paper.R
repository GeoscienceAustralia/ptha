#
# Local-scale evaluation of earthquake rates
#
library(rptha)

VARIABLE_MU = TRUE

#' Compute the PTHA18 modelled scenario rates
#'
#' Useful for regional comparisons with other studies
#'
#' @param nc_file File with the rate information (e.g.
#' ....../all_stochastic_slip_earthquake_events.nc)
#' @param segment_range the along-strike indices defining the
#' 'region-for-comparison'
#' @param num_downdip How many down-dip unit-sources on the source-zone. FIXME:
#' In future read this from a file! 
#' @param threshold_fraction_in_segment. If not null, count events where the
#' given fraction of the summed slip is in the region-for-comparison. If not
#' provided, we include all events, weighted by the fraction of their slip in
#' the region-for-comparison.
#' @param rate_of_full_length_rupture If not FALSE, then assign a weight of 1
#' to events that 'cover' the region-for-comparison along-strike, and 0 to
#' other events. This can be useful since some paleo studies report on 'full
#' length rupture' without giving magnitude specifically.
#' @param variable_mu logical. Use variable rigidity?
#' @return Look below!
get_rates_geometric_median_interpolation<-function(nc_file, segment_range, num_downdip, reference_magnitude,
    threshold_fraction_in_segment = NULL, rate_of_full_length_rupture=FALSE,
    variable_mu=VARIABLE_MU
    ){

    # Must only have 1 decimal place, to support our hard-coded geometric-mean interpolation.
    # This is identical to interpolation in log space at a halfway point between 2 values.
    # Helps deal with the fact that e.g. events with Mw=9 actually represent events in range 8.95 -- 9.05
    # Turns out we want to evaluate the rate function at '1-decimal-place' values, so some interpolation
    # is required
    stopifnot(round(reference_magnitude, 1) == reference_magnitude)

    fid = nc_open(nc_file)

    if(variable_mu){
        event_Mws = ncvar_get(fid, 'variable_mu_Mw')
        event_rate = ncvar_get(fid, 'variable_mu_rate_annual')
        event_rate_lower = ncvar_get(fid, 'variable_mu_rate_annual_lower_ci')
        event_rate_upper = ncvar_get(fid, 'variable_mu_rate_annual_upper_ci')
        event_rate_16pc = ncvar_get(fid, 'variable_mu_rate_annual_16pc')
        event_rate_84pc = ncvar_get(fid, 'variable_mu_rate_annual_84pc')
        event_rate_median = ncvar_get(fid, 'variable_mu_rate_annual_median')
    }else{
        event_Mws = ncvar_get(fid, 'Mw')
        event_rate = ncvar_get(fid, 'rate_annual')
        event_rate_lower = ncvar_get(fid, 'rate_annual_lower_ci')
        event_rate_upper = ncvar_get(fid, 'rate_annual_upper_ci')
        event_rate_16pc = ncvar_get(fid, 'rate_annual_16pc')
        event_rate_84pc = ncvar_get(fid, 'rate_annual_84pc')
        event_rate_median = ncvar_get(fid, 'rate_annual_median')
    }

    event_index_string = ncvar_get(fid, 'event_index_string')
    event_slip_string = ncvar_get(fid, 'event_slip_string')

    nc_close(fid)

    unit_source_statistics_file = Sys.glob(paste0(dirname(nc_file), '/unit_source_statistics*.nc'))
    uss = read_table_from_netcdf(unit_source_statistics_file)

    # Useful to get lon/lat range of the included unit-sources
    coordinate_corner<-function(central_lon, central_lat, uss_width, uss_length, uss_strike, uss_dip, backward){
        # Move to the top of the rupture
        p1 = destPoint(c(central_lon, central_lat), b=uss_strike-90, d=uss_width/2*cos(uss_dip/180*pi), f=0)
        if(backward){
            p2 = destPoint(p1 , b=uss_strike-180, d=uss_length/2, f=0)
        }else{
            p2 = destPoint(p1 , b=uss_strike, d=uss_length/2, f=0)
        }
        return(c(p2[1], p2[2]))
    }

    i = which(uss$alongstrike_number == segment_range[1] & uss$downdip_number == 1)
    lower_coord = coordinate_corner(uss$lon_c[i], uss$lat_c[i], uss$width[i]*1000, uss$length[i]*1000, 
                                    uss$strike[i], uss$dip[i], backward=TRUE)
    i = which(uss$alongstrike_number == segment_range[2] & uss$downdip_number == 1)
    upper_coord = coordinate_corner(uss$lon_c[i], uss$lat_c[i], uss$width[i]*1000, uss$length[i]*1000, 
                                    uss$strike[i], uss$dip[i], backward=FALSE)

    fraction_in_segment_range = mapply(f<-function(event_index_string, event_slip_string){

            inds = as.numeric(strsplit(event_index_string, '-')[[1]])
            slips = as.numeric(strsplit(event_slip_string, '_')[[1]])

            alongstrike_inds = ceiling(inds/num_downdip)

            if(!rate_of_full_length_rupture){
                k1 = sum(slips)
                k2 = sum(slips * (alongstrike_inds >= segment_range[1]) * (alongstrike_inds <= segment_range[2]))

                # Alternative approach -- fraction of area inside
                #return(mean(alongstrike_inds >= segment_range[1] & alongstrike_inds <= segment_range[2]))
                 
                return(k2 / k1)
            }else{
                # Just report 1 for full length and 0 otherwise
                k1 = as.numeric(
                    min(alongstrike_inds) <= segment_range[1] & 
                    max(alongstrike_inds) >= segment_range[2])
                return(k1)
            }
        },
        event_index_string, event_slip_string
    )


    # There are different ways of including events. We could weight every scenario rate
    # by its fraction in the segment. Or we could only include events where the fraction
    # exceeds some threshold (e.g. 0.8, or something else).
    if(is.null(threshold_fraction_in_segment)){
        fisr = fraction_in_segment_range
    }else{
        fisr = (fraction_in_segment_range > threshold_fraction_in_segment)
    }

    # NOTE: For fixed rigidity, our events with (say) Mw = 9.0 really represent
    # events with 8.95 <= Mw < 9.05
    # With variable rigidity, the magnitudes are already "smeared" both inside and outside
    # the bin. Thus, in a situation where the 'magnitude-perurbation due to variable rigidity'
    # was symmetric, then even if the scenario rates did not change, 

    r_mean  = sqrt( 
        sum((event_Mws >= reference_magnitude - 0.05)*fisr*event_rate) * 
        sum((event_Mws >= reference_magnitude + 0.05)*fisr*event_rate) )
    r_lower = sqrt(
        sum((event_Mws >= reference_magnitude - 0.05)*fisr*event_rate_lower)*
        sum((event_Mws >= reference_magnitude + 0.05)*fisr*event_rate_lower))
    r_upper = sqrt(
        sum((event_Mws >= reference_magnitude - 0.05)*fisr*event_rate_upper)*
        sum((event_Mws >= reference_magnitude + 0.05)*fisr*event_rate_upper))
    r_16pc  = sqrt(
        sum((event_Mws >= reference_magnitude - 0.05)*fisr*event_rate_16pc)*
        sum((event_Mws >= reference_magnitude + 0.05)*fisr*event_rate_16pc))
    r_84pc  = sqrt(
        sum((event_Mws >= reference_magnitude - 0.05)*fisr*event_rate_84pc)*
        sum((event_Mws >= reference_magnitude + 0.05)*fisr*event_rate_84pc))
    r_50pc  = sqrt(
        sum((event_Mws >= reference_magnitude - 0.05)*fisr*event_rate_median)*
        sum((event_Mws >= reference_magnitude + 0.05)*fisr*event_rate_median))

    #out = c(r_lower, r_16pc, r_50pc, r_mean, r_84pc, r_upper)
    output = data.frame(rp_025 = 1/r_lower, rp_mean=1/r_mean, rp_975=1/r_upper, 
                        segment_lower=segment_range[1], segment_upper=segment_range[2],
                        lower_lon = lower_coord[1], lower_lat = lower_coord[2], 
                        upper_lon = upper_coord[1], upper_lat = upper_coord[2])

    return(output)
}


get_results<-function(){

    result = list()

    #
    # Cascadia
    #
     
    segment_range = c(1, 22) # Include everything!
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/cascadia/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_cascadia.nc'
    num_downdip = 3
    reference_magnitude = 9.0 # Rong et al 2014 estimate this
    result$cascadia = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    #
    # Alaska
    #

    # Prince-William to Shumagin inclusive.
    segment_range = c(1, 29) # 
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_alaskaaleutians.nc'
    num_downdip = 4
    reference_magnitude = 9.0
    result$alaska_pw2shumagin = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    # Rate of events in the eastern region for comparison with Wesson et al. (2008)
    # The extent is based on their Table 1, including Kodiak and Prince-William-Sound segments
    # Consider extending it to include their "high-coupled" region, i.e. also with Semidi.
    # Note that the PTHA18 rates in this region are smaller than in Davies17 -- because the latter
    # used GEM's (high) coupling values, whereas PTHA18 permits low coupling.
    segment_range = c(1, 18) # 
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_alaskaaleutians.nc'
    num_downdip = 4
    reference_magnitude = 9.0
    result$alaska_wesson = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    # Butler (2016) assume length of 800km, so we slightly shortern the region
    # CONSIDER 2-17 because 17 covers the west end of Kodiak Island, which is in their data.
    # OR 3-18.
    #segment_range = c(1, 16) # 
    segment_range = c(3,18)
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_alaskaaleutians.nc'
    num_downdip = 4
    reference_magnitude = 9.0
    result$alaska_butler = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    #
    # Tohoku area
    #
    # Kagan and Jackson 2013 define this as 35-40 N, 140-146E
    segment_range = c(48, 59) # This covers kj2013, and the length (600km) is consistent with Butler2016
    num_downdip = 4
    reference_magnitude = 9.0 
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_kurilsjapan.nc'
    result$tohoku_kagan = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    # For Butler2016, same region applies
    result$tohoku_butler = result$tohoku_kagan

    #
    # Kamchatka area 
    #
    segment_range = c(8, 17) # Length consistent with Butler (500km), location based on Johnson's (1994) inversion of 1952 event
    num_downdip = 4
    reference_magnitude = 9.0
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_kurilsjapan.nc'
    result$kamchatka_butler = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    #
    # Chile -- matching 1960 earthquake
    #
    # Note Butler et al (2016) "The mean interval between the 4 most likely Mw 9 events is 667 years"
    # (but they then go on to say that because there are many more Mw 8 than Mw 9, likely some of these
    # are Mw 8).
    #
    segment_range = c(35, 53) # Extent of the 1960 event from Fuji/Satake. Consistent with Butler length of 950km
    num_downdip = 4
    reference_magnitude = 9.0
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
    result$chile_butler = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    # Mw 8.6 for comparison with Moernaut
    segment_range = c(35, 53) # Extent of the 1960 event from Fuji/Satake
    num_downdip = 4
    reference_magnitude = 8.6 
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
    result$chile_monernaut = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)


    # Mw 8.5 for comparison with Burbidge
    segment_range = c(35, 53) # Extent of the 1960 event from Fuji/Satake
    num_downdip = 4
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
    reference_magnitude = 8.5 
    result$chile_burbidge85 = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    # Mw 8 for comparison
    segment_range = c(35, 53) # Extent of the 1960 event from Fuji/Satake
    num_downdip = 4
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
    reference_magnitude = 8. 
    result$chile_burbidge8 = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    #
    # Sumatra-Andaman region
    #
    segment_range = c(69, 98) # Similar to 2004 rupture. Length consistent with Butler2016
    num_downdip = 4
    reference_magnitude = 9.0
    nc_file = '../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_sunda2.nc'
    result$sumatra_butler = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)


    #
    # Nankai
    #
    segment_range = c(1,14) #c(1,18) # Offshore of japan
    num_downdip = 3
    reference_magnitude = 8.6
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/ryuku/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_ryuku.nc'
    result$ryuku_burbidge86 = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)


    segment_range = c(1,14) #c(1,18) # Offshore of japan
    num_downdip = 3
    reference_magnitude = 8.
    nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/ryuku/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_ryuku.nc'
    result$ryuku_burbidge8 = get_rates_geometric_median_interpolation(nc_file, segment_range, num_downdip, reference_magnitude)

    return(result)
}

VARIABLE_MU=TRUE
result_variable_mu = get_results()


VARIABLE_MU=FALSE
result_constant_mu = get_results()

# Convenience function to convert to text that can go in a table
to_text_table<-function(result_site){
    
    coordinates_text = paste0(
        '(', round(result_site$lower_lon,2), ',', round(result_site$lower_lat,2), ')',
        ' to ',
        '(', round(result_site$upper_lon,2), ',', round(result_site$upper_lat,2), ')'
        )

    rp_text = paste0( round(result_site$rp_mean), ' (', round(result_site$rp_975), ',', round(result_site$rp_025), ')')

    return(c(coordinates_text, rp_text))
}

save.image('../PLOT_DATA/earthquake_rate_comparisons_paper.RData')


lapply(result_variable_mu, to_text_table)
lapply(result_constant_mu, to_text_table)

#> lapply(result_variable_mu, to_text_table)
#$cascadia
#[1] "(-125.03,40.64) to (-127.94,49.61)" "2174 (672,Inf)"                    
#
#$alaska_pw2shumagin
#[1] "(-144.49,59.21) to (-162.3,53.21)" "707 (314,5074)"                   
#
#$alaska_wesson
#[1] "(-144.49,59.21) to (-153.97,55.26)" "1230 (554,7499)"                   
#
#$alaska_butler
#[1] "(-145.3,59.15) to (-153.97,55.26)" "1319 (594,8071)"                  
#
#$tohoku_kagan
#[1] "(144.28,40.09) to (142.28,35.02)" "552 (315,1501)"                  
#
#$tohoku_butler
#[1] "(144.28,40.09) to (142.28,35.02)" "552 (315,1501)"                  
#
#$kamchatka_butler
#[1] "(162,52.48) to (157.23,48.8)" "717 (409,2205)"              
#
#$chile_butler
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "532 (313,4593)"                    
#
#$chile_monernaut
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "206 (122,678)"                     
#
#$chile_burbidge85
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "166 (98,505)"                      
#
#$chile_burbidge8
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "59 (35,142)"                       
#
#$sumatra_butler
#[1] "(96.12,1.57) to (92.68,14.45)" "728 (331,2910)"               
#
#$ryuku_burbidge86
#[1] "(138.19,34.01) to (132.43,30.81)" "1085 (380,164005)"               
#
#$ryuku_burbidge8
#[1] "(138.19,34.01) to (132.43,30.81)" "163 (87,577)"                    
#
#> lapply(result_constant_mu, to_text_table)
#$cascadia
#[1] "(-125.03,40.64) to (-127.94,49.61)" "1586 (542,Inf)"                    
#
#$alaska_pw2shumagin
#[1] "(-144.49,59.21) to (-162.3,53.21)" "770 (339,6916)"                   
#
#$alaska_wesson
#[1] "(-144.49,59.21) to (-153.97,55.26)" "1399 (628,11244)"                  
#
#$alaska_butler
#[1] "(-145.3,59.15) to (-153.97,55.26)" "1494 (670,11978)"                 
#
#$tohoku_kagan
#[1] "(144.28,40.09) to (142.28,35.02)" "649 (366,1730)"                  
#
#$tohoku_butler
#[1] "(144.28,40.09) to (142.28,35.02)" "649 (366,1730)"                  
#
#$kamchatka_butler
#[1] "(162,52.48) to (157.23,48.8)" "832 (465,2643)"              
#
#$chile_butler
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "597 (350,5889)"                    
#
#$chile_monernaut
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "216 (131,677)"                     
#
#$chile_burbidge85
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "172 (103,501)"                     
#
#$chile_burbidge8
#[1] "(-75.86,-45.62) to (-74.57,-37.58)" "60 (35,138)"                       
#
#$sumatra_butler
#[1] "(96.12,1.57) to (92.68,14.45)" "758 (353,5054)"               
#
#$ryuku_burbidge86
#[1] "(138.19,34.01) to (132.43,30.81)" "651 (268,Inf)"                   
#
#$ryuku_burbidge8
#[1] "(138.19,34.01) to (132.43,30.81)" "130 (66,481)"                    
#
