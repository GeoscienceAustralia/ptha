source('/g/data/w85/tsunami/CODE/gadi/ptha_mm/propagation/SWALS/plot.R')

# For each historic event, provide PNG name and multidomain_directory
plot_list = list(
    'tohoku_2011' = list(
        png='tohoku_2011.png', 
        md_dir='../../OUTPUTS/tohoku2011-full-ambient_sea_level_0/RUN_20241120_163301922',
        zlim=c(0,0.5)
    ),
    'solomon_2007' = list(
        png='solomon_2007.png', 
        md_dir='../../OUTPUTS/solomon2007_1_19-full-ambient_sea_level_0/RUN_20241120_163301863',
        zlim=c(0,0.5)
    )
)

# Color ramp for image
my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'), bias=2)(1000)


source('cities.R')
plot_cities <-function(){
    points(cities$lon, cities$lat, pch=19)
    text(cities$lon, cities$lat, cities$name, pos=1, cex=1.2)
}


# Large scale plots
xlim=c(110, 270)
ylim=c(-79, 60)
for(i in 1:length(plot_list)){
    png(plot_list[[i]]$png, width=7.5, height=6, units='in', res=300)
    multidomain_image(
        plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=xlim, ylim=ylim, zlim=plot_list[[i]]$zlim, 
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        NA_if_max_flux_is_zero=TRUE, # Use this for displaying a time-varying source to prevent max-stage in dry subsiding regions from being plotted.
        use_fields=TRUE
    )
    dev.off()
}

# Regional zoom plot
xlim=c(150, 154)
ylim=c(-26, -22)
for(i in 1:length(plot_list)){
    png(gsub('.png', '_region.png', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(
        plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=xlim, ylim=ylim, zlim=plot_list[[i]]$zlim, 
        asp=cos(-30/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
    )
    plot_cities()
    dev.off()
}

# Gladstone region plot
xlim=c(151.1, 151.5)
ylim=c(-24, -23.67)
for(i in 1:length(plot_list)){
    png(gsub('.png', '_gladstone.png', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(
        plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=xlim, ylim=ylim, zlim=plot_list[[i]]$zlim, 
        asp=cos(-32/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
    )
    plot_cities()
    dev.off()
}

# Yeppoon region plot
xlim=c(150.7, 151)
ylim=c(-23.33, -23.0)
for(i in 1:length(plot_list)){
    png(gsub('.png', '_yeppoon.png', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(
        plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=xlim, ylim=ylim, zlim=plot_list[[i]]$zlim, 
        asp=cos(-33/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
    )
    plot_cities()
    dev.off()
}

# Agnes waters region plot
xlim=c(151.8, 152.0)
ylim=c(-24.25, -24.13)
for(i in 1:length(plot_list)){
    png(gsub('.png', '_agnes_waters.png', plot_list[[i]]$png), width=7.5, height=6, units='in', res=300)
    multidomain_image(
        plot_list[[i]]$md_dir,
        variable='max_stage', time_index=NA, xlim=xlim, ylim=ylim, zlim=plot_list[[i]]$zlim,
        asp=cos(-28/180*pi),
        cols=my_col,
        clip_to_zlim=TRUE,
        NA_if_stage_not_above_elev=TRUE, # For consistent treatment of dry areas
        use_fields=TRUE
    )
    plot_cities()
    dev.off()
}
