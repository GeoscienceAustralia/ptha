source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')

# For each model, provide PNG name and multidomain_directory
plot_list = list(
    # Sumatra 2004
    'Sumatra2004_greater_perth_revised' = list(
        png='Sumatra2004_max_stage_greater_perth_revised2023.png', 
        md_dir='../../../greater_perth_revised2023/swals/OUTPUTS/Sumatra2004_FujiiSatake2007_timevarying-full-ambient_sea_level_0.0/RUN_20230831_172007573/',
        zlim=c(0,1)
        ),
    'Sumatra2004_bunbury_busselton' = list(
        png='Sumatra2004_max_stage_BunburyBusselton.png', 
        md_dir='../../../bunbury_busselton/swals/OUTPUTS/Fuji_andaman2004_24hrs_OpenFloodgate_NewVasseDrain_domain301122-full-ambient_sea_level_0.0/RUN_20240427_091057133/',
        zlim=c(0,1)
        ),
    'Sumatra2004_midwest' = list(
        png='Sumatra2004_max_stage_midwest.png', 
        md_dir='../../../midwest/swals/OUTPUTS/Sumatra2004_FujiiSatake2007_timevarying-full-ambient_sea_level_0.0/RUN_20231128_100150493/',
        zlim=c(0,1)
        ),
    # Sumatra 2005
    'Sumatra2005_greater_perth_revised' = list(
        png='Sumatra2005_max_stage_greater_perth_revised2023.png', 
        md_dir='../../../greater_perth_revised2023/swals/OUTPUTS/Sumatra2005_FujiiEtAl2020-full-ambient_sea_level_0.0/RUN_20230831_172051332/',
        zlim=c(0,0.3)
        ),
    'Sumatra2005_bunbury_busselton' = list(
        png='Sumatra2005_max_stage_BunburyBusselton.png', 
        md_dir='../../../bunbury_busselton/swals/OUTPUTS/Fuji_sumatra2005_newVasseData_24hrs_domain301122-full-ambient_sea_level_0.0/RUN_20240427_091106283/',
        zlim=c(0,0.3)
        ),
    'Sumatra2005_midwest' = list(
        png='Sumatra2005_max_stage_midwest.png', 
        md_dir='../../../midwest/swals/OUTPUTS/Sumatra2005_FujiiEtAl2020-full-ambient_sea_level_0.0/RUN_20231127_151302826/',
        zlim=c(0,0.3)
        )
    )

# Color ramp for image
my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'), bias=2)(1000)

# Large scale plots
for(i in 1:length(plot_list)){
    png(plot_list[[i]]$png, width=7.5, height=6, units='in', res=300)
    multidomain_image(plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=c(16, 156), ylim=c(-75, 33), zlim=plot_list[[i]]$zlim, 
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        NA_if_max_flux_is_zero=TRUE, # Use this for displaying a time-varying source to prevent max-stage in dry subsiding regions from being plotted.
        use_fields=TRUE
        )
    dev.off()
}

# Regional zoom plot
for(i in 1:length(plot_list)){
    png(paste0('Zoom_', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=c(109, 116), ylim=c(-34, -27), zlim=plot_list[[i]]$zlim, 
        asp=cos(-30/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
        )
    dev.off()
}

# Perth region plot
for(i in 1:length(plot_list)){
    png(paste0('Zoom_Perth_', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=c(115.4, 116), ylim=c(-32.32, -31.92), zlim=plot_list[[i]]$zlim, 
        asp=cos(-32/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
        )
    dev.off()
}

# Bunbury region plot
for(i in 1:length(plot_list)){
    png(paste0('Zoom_BunburyBusselton_', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=c(115, 116), ylim=c(-33.66, -32.9), zlim=plot_list[[i]]$zlim, 
        asp=cos(-33/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
        )
    dev.off()
}

# Geraldton region plot
for(i in 1:length(plot_list)){
    png(paste0('Zoom_Geraldton_', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=c(113.4, 114.9), ylim=c(-29.2, -27.7), zlim=plot_list[[i]]$zlim, 
        asp=cos(-28/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
        )
    dev.off()
}

