colorless = rgb(255, 255, 255, alpha=0, maxColorValue=255)
pal0 = cpt("ukmo_wow_humidity", n=255)[1:220]
pal1 = colorRampPalette(pal0, bias=1.8)(1024)
tmp = col2rgb(pal1)
pal1_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=160, maxColorValue=255)
# Another alternative
palDS = c('steelblue1', 'blue1', 'blue4')
palDS1 = colorRampPalette(palDS, bias=3)(1024)
tmp = col2rgb(palDS1)
palDS1_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=150, maxColorValue=255)
# Yet another alternative
palDS2 = cpt("arendal_temperature", n=400)
tmp = col2rgb(palDS2)
palDepthStage_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
palDepthStage_trans[1] = colorless
# Yet another alternative
palDS3 = cpt("cb_seq_YlOrRd_05", n=400)
tmp = col2rgb(palDS3)
palDS3_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
# And another
palSpeed = cpt("h5_jet", n=400)
tmp = col2rgb(palSpeed)
palSpeed_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
palSpeed_trans[1] = colorless
# And another
palFlux = colorRampPalette(cpt("oc_sst", n=400), bias=1.8)(400) # Concentrate variation around low values
tmp = col2rgb(palFlux)
palFlux_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
palFlux_trans[1] = colorless

# Palette for difference in max stage, zeroed on +0.0m, centre is colourless
palDiff = colorRampPalette(c('darkblue', 'blue', 'aquamarine', 'white', 'orange', 'red', 'darkred'))(1024)
palDiff[512] = colorless


# Colours for speed that will split it into 4 hazard categories using zlim=c(0, 6)
palSpeedHazard = rep(c('beige', 'orange', 'red', 'black'), each=100)
tmp = col2rgb(palSpeedHazard)
palSpeedHazard_trans = rgb(tmp[1,], tmp[2,], tmp[3,], alpha=200, maxColorValue=255)
# make the first colour transparent
palSpeedHazard_trans[1] = colorless

#
# Classes to colour the inundation rates
#
inundation_classmat = matrix(
    # Lower, upper, integerID
    c(1    ,   999999999999999,     1,
      1/100  , 1,       2, 
      1/500  , 1/100  , 3, 
      1/2500 , 1/500  , 4, 
      1/10000, 1/2500 , 5, 
      -1     , 1/10000, 6), 
    ncol=3, byrow=TRUE)
inundation_classmat_labels = c('>1', '1 - 1/100', '1/500 - 1/100', '1/2500 - 1/500', '1/10000 - 1/2500', '< 1/10000')
# Colours with transparency
inundation_classmat_col = c('cornflowerblue','red', 'orange', 'yellow', 'grey', NA)
icct = col2rgb(inundation_classmat_col, alpha=TRUE)
# the blue should be slightly transparent, and the grey
icct[4,1] = 0.7*icct[4,1]
icct[4,5] = 0.7*icct[4,5]
inundation_classmat_col = rgb(icct[1,], icct[2,], icct[3,], icct[4,], maxColorValue=255)
