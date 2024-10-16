# Append some output products that were manually created by merging the results of several models.
# This script is separate to the regular "collate_outputs.R" script because the
# files copied here were created using scripts that refer to outputs that the
# latter script creates. So it must be separate.

system('cp -r ../combined_raster_output combined_models ')
