library(cptcity)
library(rptha)

# Get the SWALS plot codes.
file_nci = '/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R'
file_home = '~/Code_Experiments/fortran/Structured_shallow_water/plot.R'
source(ifelse(file.exists(file_nci), file_nci, file_home))

# Get a 'zero-contour' [actually this is from the openstreetmap coastline dataset]
zc = readOGR('../elevation/Tonga_coast/Tonga_coast_nearlon180.shp', 'Tonga_coast_nearlon180')

all_multidomain_dirs = Sys.glob('OUTPUTS/*/RUN*')
#
# Make some plots of max-stage and elevation at various spatial scales
#
plot_max_stage_and_elevation<-function(md_dir){


    # Add in domain boundaries
    md_bbox = get_domain_interior_bbox_in_multidomain(md_dir)
    add_bboxes<-function(LWD=1, col='grey'){

        for(i in 1:length(md_bbox$merged_domain_interior_bbox)){

            dx = md_bbox$merged_domain_dx[[i]][1]
            # Avoid boxes for the 1 arcmin domains
            if(dx > 0.9*1/60) next

            LTY=c('solid', 'solid') #c('dashed', 'solid')[(dx < 0.9*1/7*1/60) + 1]
            polygon(md_bbox$merged_domain_interior_bbox[[i]], border=col, lty=LTY, lwd=LWD)
        }
    }

    # Match global domain size
    XLIM = range(md_bbox$merged_domain_interior_bbox[[1]][,1])
    YLIM = range(md_bbox$merged_domain_interior_bbox[[1]][,2])
    x_fact = 0.88

    #
    # Plot maximum stage
    #
    my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'))(1000)
    # Use a log10 transformation to stretch the max-stage colours
    stage_zlim = c(1.0e-03, 20)
    var_transform_fun<-function(x) log10(pmin(pmax(x, stage_zlim[1]), stage_zlim[2]))
    fields_axis_args = list(at=seq(log10(stage_zlim[1]),log10(stage_zlim[2])),
                            labels=round(10**(seq(log10(stage_zlim[1]),log10(stage_zlim[2]))), 3))
    out_file = paste0(md_dir, '/max_stage_full_model.png')
    png(out_file, width=8*diff(XLIM)/diff(YLIM)*x_fact, height=8, units='in', res=212.504)
    multidomain_image(md_dir, variable='max_stage', time_index=NULL,
        xlim=XLIM, ylim=YLIM, zlim=log10(stage_zlim),
        cols=my_col, clip_to_zlim=FALSE, use_fields=TRUE,
        NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE,
        var_transform_function=var_transform_fun,
        fields_axis_args=fields_axis_args)
    add_bboxes()
    dev.off()



    #
    # Plot elevation
    #

    # Colours for elevation
    #my_col = cpt("gmt_GMT_relief", n=1000)
    #ZLIM = c(-10e+03, 10e+03) # This ensures the shoreline is close to a rapid transition in the colours

    my_col = c(cpt("tp_tpsfhm", n=800), cpt("jjg_dem_c3t1", n=500))
    ZLIM = c(-8000, 5000)

    out_file = paste0(md_dir, '/elevation0_full_model.png')
    png(out_file, width=8*diff(XLIM)/diff(YLIM)*x_fact, height=8, units='in', res=212.504)
    multidomain_image(md_dir, variable='elevation0', time_index=NULL,
        xlim=XLIM, ylim=YLIM, zlim=ZLIM,
        cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE)
    add_bboxes(col='red')
    dev.off()

    # Regional domain
    XLIM = range(md_bbox$merged_domain_interior_bbox[[2]][,1])
    YLIM = range(md_bbox$merged_domain_interior_bbox[[2]][,2])
    ZLIM = c(-8000, 5000)
    x_fact = 1.2
    my_col = c(cpt("tp_tpsfhm", n=800), cpt("jjg_dem_c3t1", n=500))
    out_file = paste0(md_dir, '/elevation0_inner_model.png')
    png(out_file, width=8*diff(XLIM)/diff(YLIM)*x_fact, height=8, units='in', res=212.504)
    multidomain_image(md_dir, variable='elevation0', time_index=NULL,
        xlim=XLIM, ylim=YLIM, zlim=ZLIM,
        cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE)
    add_bboxes(col='red')
    dev.off()

    # Near Tongatapu
    XLIM = range(md_bbox$merged_domain_interior_bbox[[3]][,1])
    YLIM = range(md_bbox$merged_domain_interior_bbox[[3]][,2])
    ZLIM = c(-80, 50)
    x_fact = 1.2
    my_col = c(cpt("tp_tpsfhm", n=80), cpt("jjg_dem_c3t1", n=50))
    out_file = paste0(md_dir, '/elevation0_inner_B_model.png')
    png(out_file, width=8*diff(XLIM)/diff(YLIM)*x_fact, height=8, units='in', res=212.504)
    multidomain_image(md_dir, variable='elevation0', time_index=NULL,
        xlim=XLIM, ylim=YLIM, zlim=ZLIM,
        cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE)
    add_bboxes(col='red')
    dev.off()
}

try_plot_max_stage_and_elevation<-function(md_dir){ try(plot_max_stage_and_elevation(md_dir))}

#
# Make a movie. Assumes time-slices were stored. Requires that the ffmpeg
# module has been loaded
#
animate_stage<-function(md_dir, ZLIM=c(-10, 10)){

    library(parallel)

    start_dir = getwd()

    # Add in domain boundaries
    md_bbox = get_domain_interior_bbox_in_multidomain(md_dir)
    output_times = get_multidomain_output_times(md_dir)

    # This can be handy
    add_bboxes<-function(LWD=1, col='grey'){

        for(i in 1:length(md_bbox$merged_domain_interior_bbox)){

            dx = md_bbox$merged_domain_dx[[i]][1]
            # Avoid boxes for the 1 arcmin domains
            if(dx > 0.9*1/60) next

            LTY=c('solid', 'solid') #c('dashed', 'solid')[(dx < 0.9*1/7*1/60) + 1]
            polygon(md_bbox$merged_domain_interior_bbox[[i]], border=col, lty=LTY, lwd=LWD)
        }
    }


    #
    # Tonga wide
    #
    anim_dir = paste0(md_dir, '/animation/tongatapu_wide')
    XLIM = range(md_bbox$merged_domain_interior_bbox[[3]][,1])
    YLIM = range(md_bbox$merged_domain_interior_bbox[[3]][,2])
    COLS = (cpt('oc_sst', n=1000))
    dir.create(anim_dir, showWarnings=FALSE, recursive=TRUE)
    # This function plots the 'i'th time-slice - so we can make the plots in parallel
    parallel_fun<-function(i){
        ti = 1e+06 + i
        output_file = paste0(anim_dir, '/plot_ti_', ti, '.png')
        png(output_file, width=8, height=8*diff(YLIM)/diff(XLIM), res=300, units='in')
        multidomain_image(md_dir, variable='stage', time_index=i, asp=1/cos(mean(YLIM)/180*pi),
            xlim=XLIM, ylim=YLIM, zlim=ZLIM, cols=COLS, clip_to_zlim=TRUE,
            use_fields=TRUE, NA_if_stage_not_above_elev=FALSE)
        plot(zc, add=TRUE)
        dev.off()
    }
    mclapply(1:length(output_times), parallel_fun, mc.cores=detectCores()/2)
    # Convert to a movie and delete the pngs. Make sure ffmpeg is loaded on NCI
    setwd(anim_dir)
    system("ffmpeg -f image2 -framerate 5 -pattern_type glob -i 'plot_ti_*.png' scenario_movie.mp4")
    system("rm *.png")
    setwd(start_dir)

    #
    # Zoom to domain 4, which contains parliament house 
    #
    XLIM = range(md_bbox$merged_domain_interior_bbox[[4]][,1])
    YLIM = range(md_bbox$merged_domain_interior_bbox[[4]][,2])
    ZLIM = c(-10, 10)
    COLS = (cpt('oc_sst', n=1000))
    anim_dir = paste0(md_dir, '/animation/domain4')
    dir.create(anim_dir, showWarnings=FALSE, recursive=TRUE)
    # This function plots the 'i'th time-slice - so we can make the plots in parallel
    parallel_fun<-function(i){
        ti = 1e+06 + i
        output_file = paste0(anim_dir, '/plot_ti_', ti, '.png')
        png(output_file, width=8, height=8*diff(YLIM)/diff(XLIM), res=300, units='in')
        multidomain_image(md_dir, variable='stage', time_index=i, asp=1/cos(mean(YLIM)/180*pi),
            xlim=XLIM, ylim=YLIM, zlim=ZLIM, cols=COLS, clip_to_zlim=TRUE,
            use_fields=TRUE, NA_if_stage_not_above_elev=FALSE)
        plot(zc, add=TRUE)
        dev.off()
    }
    mclapply(1:length(output_times), parallel_fun, mc.cores=detectCores()/2)
    # Convert to a movie, and delete the pngs
    setwd(anim_dir)
    system("ffmpeg -f image2 -framerate 5 -pattern_type glob -i 'plot_ti_*.png' scenario_movie.mp4")
    system("rm *.png")
    setwd(start_dir)

}

plot_max_depth<-function(md_dir){

    all_depth_tifs = Sys.glob(paste0(md_dir, '/depth*.tif'))

    my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'))(1000)
    depth_zlim = c(1.0e-02, 10)
   
    XLIM = c(184.75, 184.9)
    YLIM = c(-21.2, -21.13) 

    output_file = paste0(md_dir, '/depth_image.png')
    ASP = 1/cos(mean(YLIM)/180*pi)
    png(output_file, width=10, height=10*diff(YLIM)/diff(XLIM) * 1/ASP * 1.5, res=300, units='in')
    plot_ext = extent(c(XLIM, YLIM))
    for(i in 3:length(all_depth_tifs)){
        r1 = raster(all_depth_tifs[i])
        r1[r1 < depth_zlim[1]] = NA
        r1[r1 > depth_zlim[2]] = depth_zlim[2]
        if(i == 3){
            plot(r1, col=my_col, asp=ASP, xlim=XLIM, ylim=YLIM, zlim=depth_zlim, ext=plot_ext)
        }else{
            image(r1, add=TRUE, zlim=depth_zlim, col=my_col)
        }
    }
    plot(zc, add=TRUE, col='slategrey', lwd=2)
    dev.off()
}
