
source('/g/data/w85/tsunami/CODE/gadi/ptha/propagation/SWALS/plot.R')

#
# Max-stage only
#
out_dir = 'multidomain_image'
dir.create(out_dir, showWarnings=FALSE)

    out_file = paste0(out_dir, '/max_stage.png')
    png(out_file, width=10, height=8, units='in', res=300)
    multidomain_image('.', variable='max_stage', time_index=NULL, 
        xlim=c(110, 300), ylim=c(-80, 65), zlim=c(0, 1), 
        cols=rev(rainbow(255)[1:200]), clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=TRUE, buffer_is_priority_domain=TRUE)
    dev.off()


    out_file = paste0(out_dir, '/max_stage_wide.png')
    png(out_file, width=17.1, height=8, units='in', res=212.504)
    multidomain_image('.', variable='max_stage', time_index=NULL, 
        xlim=c(-40, 320), ylim=c(-80, 65), zlim=c(0, 1), 
        cols=rev(rainbow(255)[1:200]), clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=TRUE, buffer_is_priority_domain=TRUE)
    dev.off()

    # As above, Greater Sydney Region
    out_file = paste0(out_dir, '/max_stage_greater_Sydney.png')
    png(out_file, width=8, height=8, units='in', res=300)
    multidomain_image('.', variable='max_stage', time_index=NULL, 
        #xlim=c(150.925, 151.532), ylim=c(-34.204, -33.406), zlim=c(0, 1), 
        xlim=c(151.0, 151.532), ylim=c(-34.100, -33.406), zlim=c(0, 1), 
        cols=rev(rainbow(255)[1:200]), clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=TRUE, buffer_is_priority_domain=TRUE, asp=1/cos(33.5/180*pi))
    dev.off()



#
# Pacific wide
#
parfun<-function(i){
    #print(i)
    out_file = paste0(out_dir, '/stage_time_index_', 1000000+i, '.png')
    png(out_file, width=10, height=8, units='in', res=300)
    multidomain_image('.', variable='stage', time_index=i, 
        xlim=c(110, 300), ylim=c(-80, 65), zlim=c(-1, 1), 
        #cols=rev(rainbow(255)[1:200]), 
        cols=grey(seq(0, 0.9, len=200)), 
	clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=TRUE, buffer_is_priority_domain=TRUE)

    if(FALSE){
	# Overlap for the chile simulaton
        par(family='serif')
        if(i > 1) text(287, -40, 'Hundreds or thousands \n (unknown) of deaths in Chile', col='darkred', adj=c(1,0), cex=1.6)
        if(i > 59) text(205, 10.5, "61 deaths \n in Hawaii", col='darkred', adj=c(0, 0), cex=1.6)
        if(i > 59) text(230, 35, "2 deaths on the \n US west coast", col='darkred', adj=c(0, 0), cex=1.6)
        if(i > 92) text(143, 30, "139 deaths \n in Japan", col='darkred', adj=c(0, 0), cex=1.6)
        if(i > 104) text(123, 10, "20 deaths in the \n Philippines", col='darkred', adj=c(0, 0), cex=1.6)
    }

    dev.off()
}
library(parallel)
mclapply(as.list(1:241), parfun, mc.cores=24)

stop('Deliberate halt')

#
# Near Australia 
#
parfun<-function(i){
    #print(i)
    out_file = paste0(out_dir, '/australia_stage_time_index_', 1000000+i, '.png')
    png(out_file, width=12, height=8, units='in', res=300)
    multidomain_image('.', variable='stage', time_index=i, 
        xlim=c(110, 190), ylim=c(-50, -5), zlim=c(-1, 1)*0.3, 
        #cols=rev(rainbow(255)[1:200]), 
        cols=grey(seq(0, 0.9, len=200)), 
	clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=TRUE, buffer_is_priority_domain=TRUE)
    dev.off()
}
library(parallel)
mclapply(as.list(1:241), parfun, mc.cores=24)

#
# SE Australia 
#
parfun<-function(i){
    #print(i)
    out_file = paste0(out_dir, '/se_australia_stage_time_index_', 1000000+i, '.png')
    png(out_file, width=12, height=8, units='in', res=300)
    multidomain_image('.', variable='stage', time_index=i, 
        xlim=c(138, 153), ylim=c(-42.5, -32.5), zlim=c(-1, 1)*0.3, 
        #cols=rev(rainbow(255)[1:200]), 
        cols=grey(seq(0, 0.9, len=200)), 
	clip_to_zlim=TRUE, use_fields=TRUE,
        NA_if_stage_not_above_elev=TRUE, buffer_is_priority_domain=TRUE)
    dev.off()
}
library(parallel)
mclapply(as.list(1:241), parfun, mc.cores=24)

