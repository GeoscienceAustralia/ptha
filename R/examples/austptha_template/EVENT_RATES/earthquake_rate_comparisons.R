#
# Local-scale evaluation of earthquake rates
#
library(rptha)

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
#' @return vector of rates with uncertainties --> c(r_2.5pc, r_16pc, r_50pc, r_mean, r_84pc, r_97.5pc)
get_rates<-function(nc_file, segment_range, num_downdip, reference_magnitude,
    threshold_fraction_in_segment = NULL, rate_of_full_length_rupture=FALSE){

    fid = nc_open(nc_file)

    event_Mws = ncvar_get(fid, 'variable_mu_Mw')

    event_rate = ncvar_get(fid, 'variable_mu_rate_annual')

    event_rate_lower = ncvar_get(fid, 'variable_mu_rate_annual_lower_ci')
    event_rate_upper = ncvar_get(fid, 'variable_mu_rate_annual_upper_ci')

    event_rate_16pc = ncvar_get(fid, 'variable_mu_rate_annual_16pc')
    event_rate_84pc = ncvar_get(fid, 'variable_mu_rate_annual_84pc')
    event_rate_median = ncvar_get(fid, 'variable_mu_rate_annual_median')

    event_index_string = ncvar_get(fid, 'event_index_string')
    event_slip_string = ncvar_get(fid, 'event_slip_string')

    nc_close(fid)


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

    r_mean  = sum((event_Mws >= reference_magnitude)*fisr*event_rate)
    r_lower = sum((event_Mws >= reference_magnitude)*fisr*event_rate_lower)
    r_upper = sum((event_Mws >= reference_magnitude)*fisr*event_rate_upper)
    r_16pc  = sum((event_Mws >= reference_magnitude)*fisr*event_rate_16pc)
    r_84pc  = sum((event_Mws >= reference_magnitude)*fisr*event_rate_84pc)
    r_50pc  = sum((event_Mws >= reference_magnitude)*fisr*event_rate_median)

    out = c(r_lower, r_16pc, r_50pc, r_mean, r_84pc, r_upper)
    return(out)
}

#
# Cascadia
#
 
segment_range = c(1, 22) # Include everything!
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/cascadia/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_cascadia.nc'
num_downdip = 3
reference_magnitude = 9.0 # Rong et al 2014 estimate this
result_cascadia = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_cascadia
#[1] 2198703.2007   63235.0438    2525.4876    2067.9503    1073.6297
#[6]     654.4533

segment_range = c(1, 22) # Include everything!
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/cascadia/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_cascadia.nc'
num_downdip = 3
reference_magnitude = 8.7 
result_cascadia = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_cascadia
#[1] 4140.9155 1692.8766  612.5835  615.8084  404.0595  278.8918

# Estimate the rate of 'full length rupture'
segment_range = c(1, 22) # Include everything!
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/cascadia/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_cascadia.nc'
num_downdip = 3
reference_magnitude = 1.0 # Just get 'full length' events
result_cascadia = get_rates(nc_file, segment_range, num_downdip, reference_magnitude, rate_of_full_length_rupture=TRUE)
#> 1/result_cascadia
#[1] 22351.3007  8117.4802  2402.7101  2124.9761  1265.9653   770.8266


#
# Alaska
#

# Rate of events in the eastern region for comparison with Wesson et al. (2007)
#
segment_range = c(1, 23) # This includes Wesson et al's Semedi segment, and east of that. 
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_alaskaaleutians.nc'
num_downdip = 4
reference_magnitude = 9.0
result_alaska = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_alaska
#[1] 3361.1557 1828.9279  997.5940  902.3920  591.5062  412.7696

#
# Butler et al. (2016). Their region might be slightly smaller, corresponding
# to 1964 event? Here we choose an extent to match Johnson's inversion. A
# little smaller than the region above.
#
segment_range = c(1, 17) 
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_alaskaaleutians.nc'
num_downdip = 4
reference_magnitude = 9.0
result_alaska = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
# [1] 4697.4216 2579.7323 1431.4851 1298.0932  854.8162  600.1659


#
# Butler et al.'s other segment range (2016)
# They call this A-A' in the paper (yellow region in Figure 1)
# Similar results to PTHA18
#
segment_range = c(18, 55) # 
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/alaskaaleutians/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_alaskaaleutians.nc'
num_downdip = 4
reference_magnitude = 9.0
result_aleutians = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_aleutians
# [1] 3127.0208 1339.5510  552.2048  506.4327  325.9272  211.4370

#
# Tohoku area
#
# Kagan and Jackson 2013 define this as 35-40 N, 140-146E
segment_range = c(48, 59)
num_downdip = 4
reference_magnitude = 9.0
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_kurilsjapan.nc'
result_tohoku = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_tohoku
#[1] 1519.8782  815.8598  580.7485  568.2235  440.8242  321.7312


#
# Kamchatka area 
#
segment_range = c(8, 17)
num_downdip = 4
reference_magnitude = 9.0
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/kurilsjapan/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_kurilsjapan.nc'
result_kamchatka = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_kamchatka
# [1] 2100.0142 1031.7235  721.9265  714.9773  553.1098  407.5389

#
# Chile -- matching 1960 earthquake
#
# Note Butler et al (2016) "The mean interval between the 4 most likely Mw 9 events is 667 years"
# (but they then go on to say that because there are many more Mw 8 than Mw 9, likely some of these
# are Mw 8).
#
segment_range = c(35, 53) # Extent of the 1960 event from Fuji/Satake
num_downdip = 4
reference_magnitude = 9.0
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
result_chile = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_chile
#[1] 1874.7230  877.3737  508.8150  531.3701  400.6907  318.1025

#Chile segment.
segment_range = c(35, 60) # Berryman et al "Central-Chile"
num_downdip = 4
reference_magnitude = 9.0
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
result_chile = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#
#> 1/result_chile
#[1] 1305.2658  603.3613  343.5816  359.7956  270.7318  213.5603

#
# Chile 2 -- chile segment -- Burbidge's estimate.
#
segment_range = c(35, 60)
num_downdip = 4
reference_magnitude = 8.5
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
result_chile = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_chile
#[1] 304.66854 175.54760 111.20742 113.29753  86.38450  68.43064

#
# Chile 3 -- chile segment -- Burbidge's estimate.
#
segment_range = c(35, 60)
num_downdip = 4
reference_magnitude = 8.0
nc_file = '/g/data/fj6/PTHA/AustPTHA_1/SOURCE_ZONES/southamerica/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_southamerica.nc'
result_chile = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#
#> 1/result_chile
#[1] 88.02668 58.13668 40.28379 40.06789 31.06519 24.54790

#
# Sumatra-Andaman region
#
segment_range = c(69, 98) # Similar to 2004 rupture
num_downdip = 4
reference_magnitude = 9.0
nc_file = '../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_sunda2.nc'
result_sumatra = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_sumatra
# [1] 2600.6321 1287.9932  783.3073  721.9693  503.7026  334.9699

segment_range = c(69, 98) # Similar to 2004 rupture
num_downdip = 4
reference_magnitude = 9.0
nc_file = '../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_sunda2.nc'
result_sumatra = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_sumatra
# [1] 2600.6321 1287.9932  783.3073  721.9693  503.7026  334.9699

# Further south, around where Patton et al.'s (2015) turbidites are
segment_range = c(48, 78)
num_downdip = 4
reference_magnitude = 9.0
nc_file = '../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_sunda2.nc'
result_sumatra = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_sumatra
#[1] 1590.7428  795.8818  515.9605  487.3140  357.0044  239.1381
reference_magnitude = 8.7
nc_file = '../SOURCE_ZONES/sunda2/TSUNAMI_EVENTS/all_stochastic_slip_earthquake_events_sunda2.nc'
result_sumatra = get_rates(nc_file, segment_range, num_downdip, reference_magnitude)
#> 1/result_sumatra
#[1] 444.9253 291.4750 212.7968 206.3841 159.5986 117.9906

