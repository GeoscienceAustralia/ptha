# Run the probabilistic inundation calculations
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4 3
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4 4
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4 5
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4 6
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4 7

# Site-specific curves at gauges
Rscript depth_vs_exrate_at_gauge.R "parliament" ptha18_tonga_MSL0_meshrefine4
Rscript depth_vs_exrate_at_gauge.R "ptha18_point_3458.3" ptha18_tonga_MSL0_meshrefine4

# Combined rasters giving exrate at given depths
Rscript raster_plots.R ptha18_tonga_MSL0_meshrefine4

# Plots
Rscript plot_depth_vs_exrate_at_parliament.R ptha18_tonga_MSL0_meshrefine4

# Offshore stage plots -- doesn't make sense to do it if MSL != 0
Rscript plot_stage_vs_exrate_at_gauge_3458.R ptha18_tonga_MSL0_meshrefine4


#
# As above with MSL = 80cm
#

# Run the probabilistic inundation calculations
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4_msl80cm 3
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4_msl80cm 4
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4_msl80cm 5
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4_msl80cm 6
Rscript probabilistic_inundation.R ptha18_tonga_MSL0_meshrefine4_msl80cm 7

# Site-specific curves at gauges
Rscript depth_vs_exrate_at_gauge.R "parliament" ptha18_tonga_MSL0_meshrefine4_msl80cm
Rscript depth_vs_exrate_at_gauge.R "ptha18_point_3458.3" ptha18_tonga_MSL0_meshrefine4_msl80cm

# Combined rasters giving exrate at given depths
Rscript raster_plots.R ptha18_tonga_MSL0_meshrefine4_msl80cm

# Plots
Rscript plot_depth_vs_exrate_at_parliament.R ptha18_tonga_MSL0_meshrefine4_msl80cm

# Offshore stage plots -- doesn't make sense to do it if MSL != 0
# Rscript plot_stage_vs_exrate_at_gauge_3458.R ptha18_tonga_MSL0_meshrefine4_msl80cm
