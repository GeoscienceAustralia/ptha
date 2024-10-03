source('../../../plot.R')


test_min_stage = function(atol=1e-6) {
    domain_number = 2
    most_recent_run = rev(Sys.glob(paste0('OUTPUTS/RUN_*')))[1]


    min_stage_raster = merge_domains_nc_grids(
        multidomain_dir = most_recent_run, 
        domain_index=domain_number,
        desired_var = 'min_stage',
        return_raster=TRUE
    )
    min_stage = as.matrix(min_stage_raster)

    # check min minimum stage is -5.987428
    check_min_stage = min(min_stage, na.rm=TRUE)
    if (abs(min(check_min_stage) - -5.987428) > atol) {
        print("FAIL")
        stop('min min_stage is not -5.987428')
    }
    print("PASS")
    # check max minimum stage is 500
    check_max_stage = max(min_stage, na.rm=TRUE)
    if (abs(check_max_stage - 500) > atol) {
        print("FAIL")
        stop('max min_stage is not 500')
    }
    print("PASS")

    # check no Inf
    if (any(is.infinite(min_stage))) {
        print("FAIL")
        stop('min_stage contains Inf')
    }
    print("PASS")

    # check no -Inf
    if (any(is.infinite(min_stage))) {
        print("FAIL")
        stop('min_stage contains -Inf')
    }
    print("PASS")

}

test_min_stage()
