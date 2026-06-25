#
# Copy information on the elevation data source
#
elevation_source_file_index_rasts = Sys.glob('../../swals/OUTPUTS/smallTestCase_kalbarri2coralbay_B_notidaladjustment-test_load_balance-ambient_sea_level_0/RUN_20260529_113644043/elevation_source_file_index_domain*.tif')
elevation_source_dir = "elevation_source_file_index"
dir.create(elevation_source_dir, showWarnings=FALSE)
file.copy(elevation_source_file_index_rasts, elevation_source_dir, overwrite=TRUE)
mydir = getwd()
setwd(elevation_source_dir)
system('gdalbuildvrt -resolution highest all_elevation_source_file_index_combined.vrt *.tif')
setwd(mydir)

# Make a legend for the elevation data source
make_source_data_legend<-function(elevation_source_dir){

    elevation_source_file = readLines('../../elevation/swals_elevation_files_in_preference_order.txt')
    esfil = data.frame(
        value = 1.0 * (1:length(elevation_source_file)), 
        quality=rep(NA, length(elevation_source_file)),
        color=rep(NA, length(elevation_source_file)),
        file = elevation_source_file) 

    # Manually populate the data quality info
    stopifnot(nrow(esfil) == 44)
    esfil$quality[1] = 'patched_manually'
    esfil$color[1] = colorRampPalette(c('purple', 'purple'))(n=1)
    esfil$quality[2:30] = 'good'
    esfil$color[2:30] = colorRampPalette(c('lightblue', 'blue'))(n=(30-2+1))
    esfil$quality[31:35] = 'regional'
    esfil$color[31:35] = colorRampPalette(c('green1', 'green4'))(n=(35-31+1))
    esfil$quality[36:44] = 'global'
    esfil$color[36:44] = colorRampPalette(c('brown2', 'brown4'))(n=(44-36+1))

    write.csv(esfil, file=paste0(elevation_source_dir, '/elevation_data_source_info.csv'), row.names=FALSE)

    # Programmatically construct the QGIS QML styling XML string
    xml_lines <- c(
      '<!DOCTYPE qgis PUBLIC "http://mrcc.com" "SYSTEM">',
      '<qgis version="3.34.0" minimumScale="0" maximumScale="1e+08" hasScaleBasedVisibilityFlag="0">',
      '  <pipe>',
      '    <rasterrenderer type="paletted" opacity="1" band="1" nodataColor="">',
      '      <colorPalette>'
    )

    # Populate the palette entries block
    for(i in 1:nrow(esfil)) {
      entry <- paste0(
        '        <paletteEntry value="', format(esfil$value[i], nsmall=1), 
        '" color="', esfil$color[i], 
        '" label="', basename(esfil$file[i]), '"/>'
      )
      xml_lines <- c(xml_lines, entry)
    }

    # Close out the XML structural tags
    xml_lines <- c(xml_lines, 
      '      </colorPalette>',
      '    </rasterrenderer>',
      '  </pipe>',
      '  <blendMode>0</blendMode>',
      '</qgis>'
    )

    # Write the lines directly to a companion file sharing the raster base name
    writeLines(xml_lines, paste0(elevation_source_dir, "/elevation_source_file_index_QGIS_raster_legend.qml"))
    return(invisible())
}
make_source_data_legend(elevation_source_dir)

