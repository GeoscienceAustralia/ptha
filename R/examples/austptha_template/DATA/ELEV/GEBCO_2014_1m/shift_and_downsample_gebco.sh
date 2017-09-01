# Set this to point to the GEBCO2014 data
export GEBCO_2014='../../../../../DATA/ELEVATION/GEBCO_2014/GEBCO_2014_2D.nc'

# Desired extends of the new elevation data
W=-39.9958333
E=320.0041667

# These parameters specify the latitude extents of the original
# file -- should be 90/-90, but since GMT uses node registration a
# half-pixel offset is added
extremeN=89.9958333
extremeS=-89.9958333

# Copy locally
cp $GEBCO_2014 GEBCO_2014_W$W-E$E-S$extremeS-N$extremeN.nc

# Change extents -- USE GMT 4 FOR THIS
grdedit -R$W/$E/$extremeS/$extremeN -S GEBCO_2014_W$W-E$E-S$extremeS-N$extremeN.nc

# Downsample from 30 arc-seconds to 1 arc-minute
gdal_translate -outsize 50% 50% -co COMPRESS=DEFLATE GEBCO_2014_W$W-E$E-S$extremeS-N$extremeN.nc GEBCO_2014_1minx1min_W$W-E$E.tif
