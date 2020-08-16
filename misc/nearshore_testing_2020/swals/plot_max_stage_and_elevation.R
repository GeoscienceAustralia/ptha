# Make a PNG with max-stage, and one with the elevation, for all models.
#
# Note these plots sometimes show linear-features or gaps at nesting boundaries
# between the outer coarse domain and the other domains. That happens because we only stored every 4'th
# grid point on the coarse domains -- so the boundaries might not match up perfectly with the
# nested domains. 
# 
# By plotting the interior bounding boxes, this is covered. Ideally we should
# just store every cell in every domain.
# 


all_multidomain_dirs = Sys.glob('OUTPUTS/*/RUN*')
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')
library(cptcity)
library(rptha)

plot_max_stage_and_elevation<-function(md_dir){

    # Add in domain boundaries
    md_bbox = get_domain_interior_bbox_in_multidomain(md_dir)
    add_bboxes<-function(LWD=1){

        for(i in 1:length(md_bbox$merged_domain_interior_bbox)){

            dx = md_bbox$merged_domain_dx[[i]][1]
            # Avoid boxes for the 1 arcmin domains
            if(dx > 0.9*1/60) next

            LTY=c('solid', 'solid') #c('dashed', 'solid')[(dx < 0.9*1/7*1/60) + 1]
            polygon(md_bbox$merged_domain_interior_bbox[[i]], border='grey', lty=LTY, lwd=LWD)
        }
    }

    #
    # Plot maximum stage
    #
    my_col = colorRampPalette(c('purple', 'blue', 'skyblue', 'green', 'yellow', 'orange', 'red', 'black'))(1000)
    # Use a log10 transformation to stretch the max-stage colours
    stage_zlim = c(1.0e-03, 1)
    var_transform_fun<-function(x) log10(pmin(pmax(x, stage_zlim[1]), stage_zlim[2]))
    fields_axis_args = list(at=seq(log10(stage_zlim[1]),log10(stage_zlim[2])),
                            labels=round(10**(seq(log10(stage_zlim[1]),log10(stage_zlim[2]))), 3))
    out_file = paste0(md_dir, '/max_stage_full_model.png')
    png(out_file, width=17.1, height=8, units='in', res=212.504)
    multidomain_image(md_dir, variable='max_stage', time_index=NULL,
        xlim=c(-40, 320), ylim=c(-79, 68), zlim=log10(stage_zlim),
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
    my_col = cpt("gmt_GMT_globe", n=1000)
    ZLIM = c(-11e+03, 6500) # This ensures the shoreline is close to a rapid transition in the colours

    out_file = paste0(md_dir, '/elevation0_full_model.png')
    png(out_file, width=17.1, height=8, units='in', res=212.504)
    multidomain_image(md_dir, variable='elevation0', time_index=NULL,
        xlim=c(-40, 320), ylim=c(-79, 68), zlim=ZLIM,
        cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE)
    add_bboxes()
    dev.off()

}

try_plot_max_stage_and_elevation<-function(md_dir){ try(plot_max_stage_and_elevation(md_dir))}

library(parallel)
MC_CORES=24
mclapply(all_multidomain_dirs, plot_max_stage_and_elevation, mc.cores=MC_CORES, mc.preschedule=FALSE)
