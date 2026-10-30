# FIXME: Add other sites.
target_sites_in_order = c(
    # "ChristmasIsland_1min_BOM",
    # "CocosIsland_BOM_1min",
    "NWONS01_02_2007_WADOT_5min",
    "EXEXM02_01_2007_WADOT_5min",
    "CNCAR02_01_2007_WADOT_5min",
    "PWLAM01_03_2007_WADOT_5min",
    "DPKBY01_03_2007a_WADOT_5min"
    # "GNGER02_02_2007a_WADOT_5min"#,
    #"Hillarys_BOM_1min_b"
)



png_name_start = 'Sumatra2007/'
dir.create(png_name_start)
png_width = 16
png_height= 6
source('plot_gauges_generic_include.R')
