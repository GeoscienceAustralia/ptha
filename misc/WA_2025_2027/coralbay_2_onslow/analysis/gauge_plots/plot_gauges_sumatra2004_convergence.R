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
png_name_start = 'Sumatra2004_convergence/'
dir.create(png_name_start)
png_width = 16
png_height= 6


plot_title_type = 'sitename'
#plot_title_type = 'lonlat'

inargs = commandArgs(trailingOnly=TRUE)

# First argument = name of RDS file containing the gauges
INPUT_RDS_COARSE = inargs[1]
INPUT_RDS_BASE = inargs[2]

# Second argument (optional) gives the YLIM scale used for all plots. If not provided
# then each gauge uses its own y-scale. Both versions are useful.
if(length(inargs) > 2){
    YLIM_RANGE = as.numeric(inargs[3])
}else{
    YLIM_RANGE = NA
}

# Third and fourth arguments (optional) give the x-axis range in hours post earthquake
if(length(inargs) > 4){
    XLIM_RANGE_HRS = as.numeric(inargs[4:5])
}else{
    XLIM_RANGE_HRS = NA
}

x_coarse = readRDS(INPUT_RDS_COARSE)
x_base = readRDS(INPUT_RDS_BASE)
xx <- list(x_base, x_coarse)
labs <- c("Baseline", "Double dx")
cols <- c("navy", "firebrick")

output_name_tag = paste0(basename(dirname(dirname(INPUT_RDS_BASE))), '_', basename(dirname(INPUT_RDS_BASE)))
if(plot_title_type == 'lonlat') output_name_tag = paste0(output_name_tag, '_gaugeLonLat_')
if(!is.na(YLIM_RANGE)) output_name_tag = paste0(output_name_tag, '_YLIM_', round(YLIM_RANGE, 4))
if(!is.na(XLIM_RANGE_HRS[1])) output_name_tag = paste0(output_name_tag, '_XLIM_', XLIM_RANGE_HRS[1], '_', XLIM_RANGE_HRS[2])

# get plot ranges from first site and first simulation
nm <- target_sites_in_order[1]
x <- xx[[1]]
model_time = x[[nm]]$model_start_julian + x[[nm]]$model_time/(3600*24)
model_stage = x[[nm]]$model_stage
obs = x[[nm]]$event_data$obs
k0 = which( (obs$juliant > min(model_time)) & (obs$juliant < max(model_time)) )
if(is.na(YLIM_RANGE)){
    YLIM = range(c(range(model_stage, na.rm=TRUE), range(obs$resid[k0], na.rm=TRUE)))
}else{
    YLIM = c(-YLIM_RANGE, YLIM_RANGE)
}
if(is.na(XLIM_RANGE_HRS)[1]){
    XLIM = range(model_time)
}else{
    XLIM = min(model_time) + XLIM_RANGE_HRS/24
}

# initialise plot
dir.create(dirname(png_name_start))
png(paste0(png_name_start, output_name_tag, '.png'), width=png_width, height=png_height, units='in', res=300)
par(mar=c(0, 3, 0, 0))
par(oma = c(4.5,1,3,0.2))
par(mfrow=c(length(target_sites_in_order), 1))

for(i in 1:length(target_sites_in_order)){

    # plot observations
    plot(XLIM, YLIM, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = XLIM, ylim = YLIM, las = 1)
    nm = target_sites_in_order[i]
    obs = x[[nm]]$event_data$obs
    obs <- obs[obs$juliant >= XLIM[1],]
    obs <- obs[obs$juliant <= XLIM[2],]
    # plot(obs$juliant, obs$resid, t='l', ylim=YLIM, xlim=XLIM, las=1, axes=FALSE)


    start_xpd = par('xpd')
    par('xpd'=NA) # Prevent clipping outside the panel
    for(j in seq_along(xx)) {
      x <- xx[[j]]

      model_time = x[[nm]]$model_start_julian + x[[nm]]$model_time/(3600*24)
      model_stage = x[[nm]]$model_stage

      # crop data to xlim
      model_stage <- model_stage[model_time >= XLIM[1]]
      model_time <- model_time[model_time >= XLIM[1]]
      model_stage <- model_stage[model_time <= XLIM[2]]
      model_time <- model_time[model_time <= XLIM[2]]

      lines(model_time, model_stage, col=cols[[j]])
    }

    # points(obs$juliant, obs$resid, t='l')
    par('xpd'=start_xpd)
    #grid(col='orange')
    grid_dy = floor(max(YLIM)*10)/10 # Nearest 10cm
    #abline(h=seq(-grid_dy, grid_dy, length=3), col='orange', lty='dashed')
    abline(h=0, col='orange', lty='dashed')
    axes_locations = x[[nm]]$model_start_julian + seq(0, 60, by=1)/24
    abline(v=axes_locations, col='orange', lty='dotted')
    axis(side=2, at=c(-grid_dy, -grid_dy/2, 0, grid_dy/2, grid_dy), line=-2, labels=FALSE)
    axis(side=2, at=c(0, grid_dy), line=-2, las=1, cex.axis=1.8)

    if(plot_title_type == 'sitename'){
        plot_title = strsplit(nm, '_')[[1]][1]
        plot_title = gsub('NWONS01', 'Beadon Creek, Onslow', plot_title)
        plot_title = gsub('EXEXM02', 'Exmouth', plot_title)
        plot_title = gsub('OnslowSalt', 'Onslow Salt', plot_title)
        plot_title = gsub('CNCAR02', 'Carnarvon', plot_title)
        plot_title = gsub('PWLAM01', 'Lower Resolution\nCape Lambert', plot_title)
        plot_title = gsub('DPKBY01', 'Lower Resolution\nKing Bay', plot_title)

        # Append the time spacing
        obs_spacing = paste0( round((obs$juliant[2] - obs$juliant[1])*24*60,1), ' min')
        plot_title = paste0(plot_title, ' (', obs_spacing, ')')
        # Remove time spacings matching (1 min) since most gauges have that, and we are short on space.
        plot_title = gsub('(1 min)', '', plot_title, fixed=TRUE) 
    }else if(plot_title_type == 'lonlat'){
        plot_title=paste0(round(x[[nm]]$model_gaugeloc[1], 3), ',', round(x[[nm]]$model_gaugeloc[2], 3))
    }


   
    text(XLIM[1], YLIM[1]*(0.5/0.7), plot_title, cex=2, pos=4)
}
axis(side=1, at=axes_locations[seq(1,61,5)], labels=seq(0, 60)[seq(1,61,5)], cex.axis=1.8)
mtext('Hours post earthquake', side=1, outer=TRUE, cex=2, line=3)
mtext('Stage (m)', side=2, outer=TRUE, cex=2, line=-1.2)
par(xpd = NA)
legend(
    "top",
    inset = c(0, -4.3),
    legend = labs,
    col = cols,
    lwd = 4,
    cex = 2,
    horiz = TRUE,
    bty = "n"
)
dev.off()
