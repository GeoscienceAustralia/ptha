# Compute the difference between two variables
# This script is used to compute the difference between two variables, and save the difference to a new raster file.
# The script will first check if the difference raster files already exist in the ./tmp directory. If not, it will compute the difference and save the result to the ./tmp directory.
library(terra)
library(parallel)


domain_range = seq(from=2, to=133)
TIF_FILES = list(sl_000=list(), sl_800=list())
TIF_FILES$sl_000$max_stage = Sys.glob(paste0('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.0/RUN_20241219_164414473/max_stage_domain_', domain_range, '.tif'))
TIF_FILES$sl_000$max_speed = Sys.glob(paste0('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.0/RUN_20241219_164414473/max_speed_domain_', domain_range, '.tif'))
TIF_FILES$sl_800$max_stage = Sys.glob(paste0('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.8/RUN_20241219_164342781/max_stage_domain_', domain_range, '.tif'))
TIF_FILES$sl_800$max_speed = Sys.glob(paste0('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.8/RUN_20241219_164342781/max_speed_domain_', domain_range, '.tif'))
ELEVATION0 = Sys.glob(paste0('../../swals/OUTPUTS/sea_level_rise/kermadec_Mw95_SLR-full-ambient_sea_level_0.0/RUN_20241219_164414473/elevation0_domain_', domain_range, '.tif'))

replace_na <- function(r, fill_raster, offset) {
    r[is.na(r)] <- fill_raster + offset
    return(r)
}

offset = list(max_stage=0.8, max_speed=0)


out_dir = './diff/'
dir.create(out_dir, showWarnings = FALSE)

for (var in c('max_stage', 'max_speed')) {
    print(var)
    n_files = length(TIF_FILES[['sl_000']][[var]])
    minimum_diff = Inf
    maximum_diff = -Inf

    # stats for computing standard deviation
    sum = 0
    sum_sq = 0
    sum_pixels = 0

    for (i in 1:n_files) {
        r1 = rast(TIF_FILES[['sl_000']][[var]][i])
        r2 = rast(TIF_FILES[['sl_800']][[var]][i])
        elev_rast = rast(ELEVATION0[i])

        to_mask = is.na(r2) | is.na(r1)
        r1 = mask(r1, to_mask)
        r2 = mask(r2, to_mask)

        if (var == 'max_speed') {
            # mask max_speed where either is NA
            tol = 1e-5
            to_mask = abs(r1) < tol | abs(r2) < tol
            r1 = mask(r1, to_mask)
            r2 = mask(r2, to_mask)
        }

        r3 = r2 - r1

        output_file = paste(out_dir, 'sl_800mm', '_', var, '_domain', i, '_diff.tif', sep='')
        writeRaster(r3, output_file, overwrite=TRUE)

        r3_matrix = as.matrix(r3)
        this_minimum_diff = min(r3_matrix, na.rm=TRUE)
        if (this_minimum_diff < minimum_diff) print(paste0(this_minimum_diff, ': ', TIF_FILES[['sl_000']][[var]][i]))
        this_maximum_diff = max(r3_matrix, na.rm=TRUE)
        if (this_maximum_diff > maximum_diff) print(paste0(this_maximum_diff, ': ', TIF_FILES[['sl_000']][[var]][i]))
        minimum_diff = min(minimum_diff, min(r3_matrix, na.rm=TRUE))
        maximum_diff = max(maximum_diff, max(r3_matrix, na.rm=TRUE))

        sum = sum + sum(sum(r3_matrix, na.rm=TRUE))
        sum_sq = sum_sq + sum(sum(r3_matrix**2, na.rm=TRUE))
        sum_pixels = sum_pixels + sum(!is.na(r3_matrix))
    }
    print(minimum_diff)
    print(maximum_diff)

    mean_ = sum / sum_pixels
    var_ = (sum_sq / sum_pixels) - (mean_**2)
    std_ = (var_)**.5

    print(paste0("Standard deviation: ", std_))
    print(paste0("Mean: ", mean_))
}
