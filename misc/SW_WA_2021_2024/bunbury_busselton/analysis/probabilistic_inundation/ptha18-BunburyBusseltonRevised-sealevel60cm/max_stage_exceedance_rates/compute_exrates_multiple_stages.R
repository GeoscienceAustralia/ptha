#
# Sum exceedance-rates for multiple max-stage values over the sunda arc
#
#
library(raster)

# We did calculations for these max-stage values
max_stages_char = paste0('_', as.character( 0.6 + c(0.001, 1:10) ), '.tif')
# And these domains -- 526 in total
domains_char = paste0('_domain_', seq(1, 526), '_')

# Store results here
output_dirs = paste0('summed_results', gsub('.tif', '', max_stages_char, fixed=TRUE))
for(i in 1:length(output_dirs)) dir.create(output_dirs[i], showWarnings=FALSE)

# tifs for all domains, all max-stages, outerrisesunda
outerrisesunda_tifs = Sys.glob('ptha18-BunburyBusseltonRevised-sealevel60cm-random_outerrisesunda-/*.tif')
# tifs for all domains, all max-stages, all segments and unsegmented model, sunda2
sunda2_tifs = Sys.glob('ptha18-BunburyBusseltonRevised-sealevel60cm-random_sunda2-/*.tif')


sum_tifs<-function(j, i){
    # Process tifs corresponding to max_stages_char[j], domains_char[i]
    # 'outerrise-sunda exrate' + 0.5 * (sum of unsegmented and segmented 'sunda2-exrates')

    output_dir = output_dirs[j]

    # Find outerrisesunda tif -- only unsegmented
    outerrise_tif_ind = which(grepl(max_stages_char[j], outerrisesunda_tifs) & 
                              grepl(domains_char[i], outerrisesunda_tifs))
    stopifnot(length(outerrise_tif_ind) == 1)

    # Find sunda2 tifs -- unsegmented + one per segment
    sunda2_tif_ind = which(grepl(max_stages_char[j], sunda2_tifs) & 
                           grepl(domains_char[i], sunda2_tifs))
    stopifnot(length(sunda2_tif_ind) == 5)

    # Process the results for sunda2 --
    # 
    sunda2_tif = sunda2_tifs[sunda2_tif_ind]
    output = raster(sunda2_tif[1])
    output[is.na(output)] = 0
    for(ii in 2:length(sunda2_tif)){
        tmp = raster(sunda2_tif[ii])
        tmp[is.na(tmp)] = 0
        # Sum here -- later multiply by 0.5
        output = output + tmp
    }
    output = 0.5 * output # 50% segmented, 50% unsegmented

    # Add on the outerrisesunda results
    outerrise_tif = outerrisesunda_tifs[outerrise_tif_ind]
    tmp = raster(outerrise_tif)
    tmp[is.na(tmp)] = 0
    output = output + tmp

    # Save result
    output_file = paste0(output_dir, '/', gsub('unsegmented', 'combined', basename(outerrise_tif)))
    writeRaster(output, file=output_file, overwrite=TRUE, options=c('COMPRESS=DEFLATE'))

    rm(output, tmp)
    gc()

    return(1)

}

# Wrap in try to avoid parallel failure
try_sum_tifs<-function(inlist){
    try(sum_tifs(inlist$j, inlist$i))
}

input_dat = expand.grid(1:length(max_stages_char), 1:length(domains_char))
inlist = vector(mode='list', length=nrow(input_dat))
for(i in 1:nrow(input_dat)){
    inlist[[i]] = list(j=input_dat[i,1], i = input_dat[i,2])
}

library(parallel)
mclapply(inlist, try_sum_tifs, mc.cores=12)
