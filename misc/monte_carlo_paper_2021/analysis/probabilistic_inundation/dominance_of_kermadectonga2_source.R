#
# At Site P, figure out how much the kermadectonga2 source contributes to the hazard, vs other sites.
#
library(rptha)

#
# Get the PTHA18 access codes, needed for functions below
#
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/ptha_access/get_PTHA_results.R'
file_home = '/media/gareth/Windows7_OS/Users/gareth/Documents/work/AustPTHA/CODE/ptha/ptha_access/get_PTHA_results.R'
ptha18 = new.env()
source(ifelse(file.exists(file_nci), file_nci, file_home), local=ptha18, chdir=TRUE)

# Get the exceedance-rate curves for all source zones
rate_at_siteP_all_sources = ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
    hazard_point_gaugeID = 3458.3,
    only_mean_rate_curve=TRUE)

rate_at_siteP_kt2 = ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
    hazard_point_gaugeID = 3458.3,
    source_name = 'kermadectonga2',
    only_mean_rate_curve=TRUE)

rate_at_siteP_outerrisekt = ptha18$get_stage_exceedance_rate_curve_at_hazard_point(
    hazard_point_gaugeID = 3458.3,
    source_name = 'outerrise_kermadectonga',
    only_mean_rate_curve=TRUE)

all_sources_rates = approx(rate_at_siteP_all_sources$stage, rate_at_siteP_all_sources$stochastic_slip_rate, xout=c(2,4))$y
kt2_rates         = approx(rate_at_siteP_kt2$stage        , rate_at_siteP_kt2$stochastic_slip_rate        , xout=c(2,4))$y
outerrise_rates   = approx(rate_at_siteP_outerrisekt$stage, rate_at_siteP_outerrisekt$stochastic_slip_rate, xout=c(2,4))$y

#
# How much does the kermadectonga2 thrust source-zone contribute?
#

#> all_sources_rates
#[1] 0.0013641317 0.0003696275
#> kt2_rates
#[1] 0.0012243411 0.0003538468
#> kt2_rates/all_sources_rates
#[1] 0.8975241 0.9573066

#
# What about the combination of the kermadectonga2 thrust and outer-rise sources?
#

#
#> (kt2_rates + outerrise_rates)/all_sources_rates
#[1] 0.9778804 1.0000000
