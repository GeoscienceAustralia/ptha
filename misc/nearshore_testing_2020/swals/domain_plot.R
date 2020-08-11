
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')
library(cptcity)
library(rptha)

#
# Elevation0 only
#
out_dir = 'multidomain_image'
dir.create(out_dir, showWarnings=FALSE)

# Zero contour
zc = readOGR('../sources/figures/zero_contour/zero_contour.shp', layer='zero_contour')

# Add in domain boundaries
md_bbox = get_domain_interior_bbox_in_multidomain('.')
add_bboxes<-function(LWD=1){

    plot(zc, add=TRUE)

    for(i in 1:length(md_bbox$merged_domain_interior_bbox)){

        dx = md_bbox$merged_domain_dx[[i]][1]
        # Avoid boxes for the 1 arcmin domains
        if(dx > 0.9*1/60) next
       
        LTY=c('dashed', 'solid')[(dx < 0.9*1/7*1/60) + 1]
        polygon(md_bbox$merged_domain_interior_bbox[[i]], border='red', lty=LTY, lwd=LWD)
    }
}


# Colours for elevation
my_col = cpt("gmt_GMT_globe", n=1000)
ZLIM = c(-11e+03, 6500) # This ensures the shoreline is close to a rapid transition in the colours

# Pacific focus
out_file = paste0(out_dir, '/elevation0_Pacific.png')
png(out_file, width=10, height=8, units='in', res=300)
multidomain_image('.', variable='elevation0', time_index=NULL, 
    xlim=c(110, 300), ylim=c(-80, 65), zlim=ZLIM, 
    cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
    NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE)
add_bboxes()
dev.off()

# Global
out_file = paste0(out_dir, '/elevation0_full_model.png')
png(out_file, width=17.1, height=8, units='in', res=212.504)
multidomain_image('.', variable='elevation0', time_index=NULL, 
    xlim=c(-40, 320), ylim=c(-80, 65), zlim=ZLIM, 
    cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
    NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE)
add_bboxes()
dev.off()

# As above, Greater Sydney Region
out_file = paste0(out_dir, '/elevation0_greater_Sydney.png')
png(out_file, width=8, height=8, units='in', res=300)
multidomain_image('.', variable='elevation0', time_index=NULL, 
    xlim=c(151.0, 151.532), ylim=c(-34.100, -33.406), zlim=ZLIM/30, 
    cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
    NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE, 
    asp=1/cos(33.5/180*pi))
add_bboxes(LWD=2)
dev.off()

# As above, Australia
out_file = paste0(out_dir, '/elevation0_Australia.png')
png(out_file, width=12, height=8, units='in', res=300)
multidomain_image('.', variable='elevation0', time_index=NULL, 
    xlim=c(110.0, 160), ylim=c(-45, -10), zlim=ZLIM/2, 
    cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
    NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE, 
    asp=1/cos(30/180*pi))
add_bboxes(LWD=2)
dev.off()

# NSW/VIC
out_file = paste0(out_dir, '/elevation0_NSWVIC.png')
png(out_file, width=12, height=8, units='in', res=300)
multidomain_image('.', variable='elevation0', time_index=NULL, 
    xlim=c(140.0, 155), ylim=c(-42, -30), zlim=ZLIM/2, 
    cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
    NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE, 
    asp=1/cos(35/180*pi))
add_bboxes(LWD=2)
dev.off()

# WA
out_file = paste0(out_dir, '/elevation0_SWWA.png')
png(out_file, width=12, height=8, units='in', res=300)
multidomain_image('.', variable='elevation0', time_index=NULL, 
    xlim=c(110.0, 120), ylim=c(-37.5, -29), zlim=ZLIM/2, 
    cols=my_col, clip_to_zlim=TRUE, use_fields=TRUE,
    NA_if_stage_not_above_elev=FALSE, buffer_is_priority_domain=FALSE, 
    asp=1/cos(32/180*pi))
add_bboxes(LWD=2)
dev.off()


add_parallel_partition<-function(){

    # All colours, no grey
    COLORS = colors()[grep('gr(a|e)y', colors(), invert = T)]
    set.seed(1234)
    my_col = sample(COLORS, size=max(md_bbox$domain_images))

    # Get an order of polygons from 'coarsest' to 'less coarse' to define plotting order
    poly_order = order(unlist(lapply(md_bbox$domain_dx, f<-function(x) x[1])), decreasing=TRUE)

    for(i in poly_order) polygon(md_bbox$domain_interior_bbox[[i]], col=my_col[md_bbox$domain_images[i]])

    plot(zc, add=TRUE)

}

out_file = paste0(out_dir, '/partition_full_model.png')
png(out_file, width=17.1, height=8, units='in', res=212.504)
plot(c(0, 1), c(0, 1), col=0,  xlim=c(-40, 320), ylim=c(-80, 65), asp=1)
add_parallel_partition()
dev.off()

out_file = paste0(out_dir, '/partition_Australia.png')
png(out_file, width=12, height=8, units='in', res=300)
plot(c(0, 1), c(0, 1), col=0, xlim=c(110.0, 160), ylim=c(-45, -10), asp=1/cos(30/180*pi))
add_parallel_partition()
dev.off()

out_file = paste0(out_dir, '/partition_NSWVIC.png')
png(out_file, width=12, height=8, units='in', res=300)
plot(c(0, 1), c(0, 1), col=0, xlim=c(140.0, 155), ylim=c(-42, -30), asp=1/cos(35/180*pi))
add_parallel_partition()
dev.off()

out_file = paste0(out_dir, '/partition_SWWA.png')
png(out_file, width=12, height=8, units='in', res=300)
plot(c(0, 1), c(0, 1), col=0, xlim=c(110.0, 120), ylim=c(-37.5, -29), asp=1/cos(35/180*pi))
add_parallel_partition()
dev.off()

