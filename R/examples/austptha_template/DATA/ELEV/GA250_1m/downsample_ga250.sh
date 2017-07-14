# Set this to point to the GA250 data file
export GA_250='../../../../../DATA/ELEVATION/GA250m/ER_Mapper_ers/ausbath_09_v4_ex_ex.ers'

# The extent is 80 deg lon x 52 deg lat, so we multiply those by 60 to get the
# desired grid size (1 arc minute)
gdal_translate -outsize 4800 3120 -co COMPRESS=DEFLATE $GA_250 GA_250_1mx1m.tif
