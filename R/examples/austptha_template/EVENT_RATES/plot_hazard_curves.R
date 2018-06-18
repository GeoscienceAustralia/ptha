source('plot_hazard_curves_utilities.R')

#
# If this code is called from Rscript, then make the plots
#
if(interactive() == FALSE){

    make_global_stage_return_period_plots()
    make_stage_exceedance_rate_convergence_plot()
    make_standard_site_exceedance_rate_plots()

}
