target_sites_in_order = c(
    #"ChristmasIsland_1min_2004",
    # "CocosIsland_BOM_1min_2004",
    "NWONS01_02_2004_WADOT_5min",
    "EXEXM02_01_2004_WADOT_5min",
    "CNCAR02_01_2004_WADOT_5min",
    "PWLAM01_03_2004_WADOT_5min",
    "DPKBY01_03_2004a_WADOT_5min"
    # "GeraldtonAlternate_WADoT_5min_2004"#,
    #"Hillarys_BOM_1min_2004"
    )
png_name_start = 'Sumatra2004/'
dir.create(png_name_start)
png_width = 16
png_height= 6
source('plot_gauges_generic_include.R')
