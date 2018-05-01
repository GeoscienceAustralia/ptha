#
# See if there are any 'general' patterns between GOF statistics and event
# properties. E.G. Maybe the 'good' events consistently have small area
# compared to the corresponding family, or low peak slip compared to the
# corresponding family, etc.
#

#' Extract statistics from the corresponding family of model scenarios
#'
#' Works on objects made in gauge_summary_statistics.R, that are lists
# (per dart buoy) containing lists (per model scenario) with summary statistics
#'
#' Returns a data.frame with the goodness of fit statistic and other columns summarising
#' the rupture.
#' 
#' @param gauge_stats object like 'stochastic_slip_stats' or 'uniform_slip_stats', etc, as created
#' by the script gauge_summary_statistics.R {in e.g. SOURCE_ZONES/sourcename/TSUNAMI_EVENTS/plots/
#' @param unit_source_statistics the unit source statistics
family_stats<-function(gauge_stats, unit_source_statistics){

    # Get time goodness-of-fit statistic for each model scenario
    gf_mat = lapply(gauge_stats, f<-function(x) lapply(x, f<-function(x) x$model_data_similarity_time))
    # Convert from list of lists to matrix
    for(i in 1:length(gf_mat)) gf_mat[[i]] = unlist(gf_mat[[i]])
    gf_mat = matrix(unlist(gf_mat), ncol=length(gf_mat))
    gf_median = apply(gf_mat, 1, median)

    # Get the peak slip for each model scenario
    if('event_slip_string' %in% names(gauge_stats[[1]][[1]]$events_with_Mw)){
        peak_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) max(as.numeric(strsplit(x$events_with_Mw$event_slip_string, '_')[[1]]))))
        # Mean slip 
        mean_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) mean(as.numeric(strsplit(x$events_with_Mw$event_slip_string, '_')[[1]]))))
    }else{
        peak_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$slip))
        mean_slip_sum = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$slip))
    }

    # Nonzero slip area 
    area = unlist(lapply(gauge_stats[[1]], f<-function(x){
            inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
            area = sum(unit_source_statistics$length[inds] * unit_source_statistics$width[inds])
            return(area)
        }))

    # Length -- we find the alongstrike range of the event, and sum the near trench lengths
    length = unlist(lapply(gauge_stats[[1]], f<-function(x){
            inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
            alongstrike_range = range(unit_source_statistics$alongstrike_number[inds])
            inds2 = which(unit_source_statistics$alongstrike_number >= alongstrike_range[1] & 
                          unit_source_statistics$alongstrike_number <= alongstrike_range[2] & 
                          unit_source_statistics$downdip_number == 1)
            length = sum(unit_source_statistics$length[inds2])
            return(length)
            }))
            
    # width -- we find the downdip range of the event, and sum the widths at the least-along-strike location
    width = unlist(lapply(gauge_stats[[1]], f<-function(x){
            inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
            downdip_range = range(unit_source_statistics$downdip_number[inds])
            inds2 = which(unit_source_statistics$downdip_number >= downdip_range[1] & 
                          unit_source_statistics$downdip_number <= downdip_range[2] & 
                          unit_source_statistics$alongstrike_number == unit_source_statistics$alongstrike_number[inds[1]])
            width = sum(unit_source_statistics$width[inds2])
            return(width)
            }))
            

    if('physical_corner_wavenumber_x' %in% names(gauge_stats[[1]][[1]]$events_with_Mw)){
        # Variable slip

        # Corner wavenumbers
        corner_wavenumber_x = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$physical_corner_wavenumber_x))
        corner_wavenumber_y = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$physical_corner_wavenumber_y))

        # Peak slip location
        peak_slip_downdip = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$peak_slip_downdip_ind))
        peak_slip_alongstrike = unlist(lapply(gauge_stats[[1]], f<-function(x) x$events_with_Mw$peak_slip_alongstrike_ind))
    }else{
        # Both uniform slip models

        corner_wavenumber_x = width*NA
        corner_wavenumber_y = width*NA

        # Peak slip location
        peak_slip_downdip = unlist(lapply(gauge_stats[[1]], f<-function(x){
                inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
                downdip_mean = mean(unit_source_statistics$downdip_number[inds])                 
                return(downdip_mean) 
        }))
        # Peak slip location
        peak_slip_alongstrike = unlist(lapply(gauge_stats[[1]], f<-function(x){
                inds = as.numeric(strsplit(x$events_with_Mw$event_index_string, '-')[[1]])
                alongstrike_mean = mean(unit_source_statistics$alongstrike_number[inds])                 
                return(alongstrike_mean) 
        }))


    }

    output = data.frame(gf = gf_median, peak_slip = peak_slip_sum, mean_slip = mean_slip_sum, area = area,
        length=length, width=width, corner_wavenumber_x = corner_wavenumber_x, corner_wavenumber_y = corner_wavenumber_y,
        peak_slip_downdip = peak_slip_downdip, peak_slip_alongstrike = peak_slip_alongstrike)

    return(output)
}

#
#
# Main script here
#
#

variable_mu = FALSE

if(variable_mu){
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*varyMu.Rdata')
}else{
    all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*[0-9].Rdata')
}

uniform_stat = vector(mode='list', length=length(all_Rdata))
stochastic_stat = vector(mode='list', length=length(all_Rdata))
variable_uniform_stat = vector(mode='list', length=length(all_Rdata))

for(i in 1:length(all_Rdata)){

    print(all_Rdata[i])
    event_env = new.env()

    # Load the R session associated with the gauge_summary_statistics.R script
    load(all_Rdata[i], envir=event_env)

    # Main computation here
    stochastic_stat[[i]] = family_stats(event_env$stochastic_slip_stats, event_env$unit_source_statistics)
    uniform_stat[[i]] = family_stats(event_env$uniform_slip_stats, event_env$unit_source_statistics)
    variable_uniform_stat[[i]] = family_stats(event_env$variable_uniform_slip_stats, event_env$unit_source_statistics)
}
event_env = new.env() # Clear the memory
names(uniform_stat) = basename(all_Rdata)
names(stochastic_stat) = basename(all_Rdata)
names(variable_uniform_stat) = basename(all_Rdata)

save.image('event_properties_and_GOF_session.Rdata')


