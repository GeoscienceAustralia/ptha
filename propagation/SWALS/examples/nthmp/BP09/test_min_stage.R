source('../../../plot.R')


test_min_stage = function(rtol=1e-2) {
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
    expected_val = -5.987428
    if (abs(min(check_min_stage) - expected_val)/abs(expected_val) > rtol) {
        print(paste0("FAIL - min min_stage is not within rtol of ", expected_val))
    }else{
        print("PASS")
    }
    # check max minimum stage is 500
    check_max_stage = max(min_stage, na.rm=TRUE)
    expected_val = 500
    if (abs(check_max_stage - expected_val)/abs(expected_val) > rtol) {
        print(paste0("FAIL: max min_stage is not within rtol of ", expected_val))
    }else{
        print("PASS")
    }

    # check no Inf
    if (any(is.infinite(min_stage))) {
        print("FAIL - min_stage contains Inf")
    }else{
        print("PASS")
    }

    # check no -Inf
    if (any(is.infinite(min_stage))) {
        print("FAIL - min_stage contains -Inf")
    }else{
        print("PASS")
    }

}

test_min_stage()
