#' Create a single initial condition raster for a given row in the PTHA18 results
#' Outputs a single raster file to the test/large_kermadectonga2_43783/ directory

ptha18 = new.env()
source(
    '/g/data/w85/tsunami/CODE/gadi/ptha_mm/ptha_access/get_PTHA_results.R',
    ptha18,
    chdir=TRUE
)

target_row = 148059
local_Mw = 9.6
source_zone = 'southamerica'
output_tif_dir = 'large_southamerica_148059/'

# Make a row index character, padded with leading zeros as required, for consistency of file names
target_row_fixed_width = substring(as.character(1e+07 + target_row), 2, 8)

dir.create(output_tif_dir, showWarnings=FALSE)
output_filename = paste0(
    output_tif_dir,
    '/',
    source_zone,
    '_row_',
    target_row_fixed_width,
    '_Mw_',
    round(local_Mw*10),
    '_HS.tif'
)

source_zone_scenario = ptha18$get_source_zone_events_data(source_zone, slip_type='stochastic')
r1 = ptha18$get_initial_condition_for_event(source_zone_scenario, event_ID=target_row)
writeRaster(r1, output_filename, options=c('COMPRESS=DEFLATE'), overwrite=TRUE)


# Eyeball the sources
plot_all_initial_conditions<-function(){

    all_ic = Sys.glob(paste0(output_tif_dir, '/*.tif'))
    pdf_file = paste0(output_tif_dir, gsub(".tif", ".pdf", basename(all_ic)))
    pdf(pdf_file, width=6, height=9)
    for(i in 1:length(all_ic)){
        r1 = raster(all_ic[i])
        plot(r1)
        title(basename(all_ic[i]))
    }
    dev.off()
}

plot_all_initial_conditions()