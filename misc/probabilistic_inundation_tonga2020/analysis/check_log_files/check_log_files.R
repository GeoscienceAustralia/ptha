# Code to check the log files, and notes about the interpretation of energy.

# 
# NOTES ON INTERPRETATION OF SWALS ENERGY CALCULATIONS
#
# As outlined below, the SWALS energy calculations can be affected by the following:
# - Boundary fluxes -- the boundary condition we used in this model is not guarenteed to
#   conserve or dissipate energy (it depends).
# - Okada initial conditions with non-zero volume -- this will make the SWALS
#   potential energy differ from the available potential energy when (MSL_LINEAR != 0)
# - Wetting and drying -- SWALS potential energy computations should be valid in this
#   case, but the classical equation {integral g/2 stage^2} will not be.
# - If we run a model at different mean-sea-levels, then it is possible to derive
#   the expected relation between the energy in the models (ignoring differences caused by the different
#   MSL). We get an equation like:
#       normalised_energy_on_rho = energy_on_rho + grav * msl_linear * (md_volume[1] - md_volume)
#   where energy_on_rho is what SWALS_reports, md_volume is the time-varying volume in the multidomain, 
#   and msl_linear is the SWALS msl_linear. This is a good way to compare the energy in models that
#   are "the same" except for different MSL_LINEAR values.
#
########
#
# THEORY
#
########
#
# Below are notes on the relation between the SWALS potential energy term, and the "available potential energy".
# The relationship is slightly complex. It is worth noting there are not such complexities for the kinetic energy,
# so that can be a more straightforward metric to analyze.
#
# In SWALS we compute potential energy as:
#    SWALS_potential_energy_on_rho = g/2 * integral{ stage^2 - elev^2 } - CONSTANT
# where:
#    CONSTANT = g/2 * integral_where_elevation_less_than_MSL_LINEAR { MSL_LINEAR^2 - elev^2 }
# is designed to make the potential energy zero when the wet-stage is equal to MSL_LINEAR everywhere.
#
# How does this compare with the classical equation for "available potential energy" which
# ignores wetting and drying?:
#    available_potential_energy_on_rho = g/2 * integral { (stage - MSL_LINEAR)^2 }
# 
# Supposing there is no wetting and drying, we can expand out the the SWALS potential energy calculation as:
#     SWALS_potential_energy_on_rho = g/2 * integral{ stage^2 - MSL_LINEAR^2}
# whereas the classical equation expands as.
#     available_potential_energy_on_rho = g/2 * integral{ stage^2 - 2*stage * MSL_LINEAR + MSL_LINEAR^2}
# If we suppose the mean_stage REALLY is MSL_LINEAR, i.e.
#    integral{stage} == integral{MSL_LINEAR}
# then 
#    integral{ -2*stage*MSL_LINEAR } = integral{ -2 * MSL_LINEAR^2 }
# so we can rearrange the classical equation as 
#     available_potential_energy_on_rho = g/2 * integral{ stage^2 - MSL_LINEAR^2} = SWALS_potential_energy_on_rho
# In summary, the SWALS calculation is the same as the classical available potential energy if there
# is no wetting and drying, and the mean-wet-stage == MSL_LINEAR.
#
# In general that might not be the case, however [e.g. if we initialise the model with an Okada initial
# condition that has non-zero volume, the mean-wet-stage will no longer be as expected.]. In that case,
# if there is no wetting and drying, then [KEY EQUATION]:
#     available_potential_energy_on_rho = SWALS_potential_energy_on_rho + g/2 * 2 * MSL_LINEAR*integral{ (MSL_LINEAR - stage) }
# Notice that the time-variation of the latter integral is related to the
# volume of water in the domain, which is something our volume conservation
# analysis keeps track of. Also, the latter term vanishes if MSL_LINEAR=0
#
# So the SWALS_potential_energy_on_rho may be greater than or less than the available_potential_energy_on_rho. We should still
# have conservation of the SWALS potential energy [ which is the actual potential energy plus some constant ], or
# dissipation in models with friction, but it does not correspond to the available potential energy (unless MSL_LINEAR=0).
# However, in the latter case, if we understand the change in volume of the domain over time, and know MSL_LINEAR, then we
# can account for those factors [as done in the function below]
#
# There is also the issue that the boundaries are open, with radiation conditions that nudge the stage
# toward the MSL_LINEAR. This means that the energy might not always decrease at later times.
# 
# Summary: Wetting and drying will invalidate the classical
# available_potential_energy calculation, but not SWALS calculation. However,
# if the mean stage is not MSL_LINEAR, then the SWALS calculation will not
# correspond to the available potential energy. Boundary fluxes make all
# reasoning about the energy balance difficult.
# 

get_time_energy_stagerange_from_log<-function(log_file){

    # Get the MSL_linear from the log file name
    split_name = strsplit(log_file, '-')
    msl_linear = as.numeric(gsub("ambientsealevel_","", split_name[[1]][3], fixed=TRUE))

    x = readLines(log_file)
    # Find log-lines just above the reported total energy
    k = grep('Global energy-total', x)
    energy_on_rho = as.numeric(x[k+1])
    kinetic_energy_on_rho = as.numeric(x[k-1])

    # In these log files, the desired time (secs) happens to be in index "k-33"
    # (detail of the sequence of print-outs, likely to change in future)
    if(all(substring(x[k-34], 1, 5) == 'Time:')){
        time = as.numeric(x[k-33])
    }else{
        stop('Did not find "Time:" in index (k-34) -- has the log file format changed?')
    }

    # Also store max-stage
    max_stage_data = as.numeric(x[k-9])
    min_stage_data = as.numeric(x[k-8])

    # Volume conservation error
    mass_error = as.numeric(substring(x[k+7], 29, nchar(x[k+7])))
    volume_change = as.numeric(substring(x[k+5], 29, nchar(x[k+5])))
    # volume_change + boundary_flux_integral = "unexplained change" (mass error)

    # Boundary_flux integral.
    # Beware sometimes ifort does strange formatting of the boundary_flux_integral line,
    # it should be like:
    #    boundary flux integral:  -1.334065992544E-99
    # but we can end up missing the "E" when the magnitude is < 1e-100:
    #     boundary flux integral:  -1.255180641882-130
    # I have since increased the digits for the exponent in the format specifiers. But a
    # work-around follows from the relation 
    #    "mass_error = volume_change + boundary_flux_integral"
    # and the fact that we don't hit the formatting issue for those runs.
    boundary_flux_integral = as.numeric(substring(x[k+6], 29, nchar(x[k+6])))
    # Here we work around the aforementioned format issue in some log files.
    alternate_boundary_flux_integral = mass_error - volume_change
    b = which(is.na(boundary_flux_integral))
    if(length(b) > 0) boundary_flux_integral[b] = alternate_boundary_flux_integral[b]

    # Volume
    md_volume =  as.numeric(substring(x[k+4], 29, nchar(x[k+4])))

    # This quantity makes it easier to compare the energy results from both models,
    # by transforming them to something closer to the available potential energy in
    # the no-wetting-drying case
    grav = 9.8
    normalised_energy_on_rho = energy_on_rho + grav * msl_linear * (md_volume[1] - md_volume)

    return(list(time=time, 
                energy_on_rho=energy_on_rho, 
                kinetic_energy_on_rho=kinetic_energy_on_rho,
                normalised_energy_on_rho = normalised_energy_on_rho,
                max_stage = max_stage_data, 
                min_stage=min_stage_data,
                log_file=log_file, 
                mass_error = mass_error, 
                boundary_flux_integral = boundary_flux_integral,
                md_volume = md_volume, 
                msl_linear=msl_linear))
}

# Logical checks on the result of get_time_energy_stagerange_from_log
check_log<-function(log_data){

    n = length(log_data$time)

    # Ad-hoc indicator for the time when boundary-fluxes are zero [i.e.  before
    # the wave reaches a model boundary]. After the model starts having
    # boundary fluxes, we can't guarentee anything about energy conservation.
    flux_threshold = 1
    before_bc = min(which(abs(log_data$boundary_flux_integral) > flux_threshold)) - 1
    if(!is.finite(before_bc)) before_bc = n

    before_bc_limit = before_bc

    if(before_bc <= 1) before_bc = 3

    water_density = 1024

    data.frame(
        n = n,
        before_bc = before_bc_limit, # Number of time-steps before the stage-perturbation reaches the boundary, with min=2
        run_finished = (log_data$time[n] == 18000),
        stage_maxima = max(log_data$max_stage),
        stage_minima = min(log_data$min_stage),
        initial_volume = log_data$md_volume[1],
        energy_start = log_data$energy_on_rho[1] * water_density,
        energy_end = log_data$energy_on_rho[n] * water_density,
        energy_max = max(log_data$energy_on_rho) * water_density,
        energy_max_while_closed = max(log_data$energy_on_rho[1:before_bc]) * water_density,
        # Th change in energy is hard to interpret once boundary fluxes can
        # happen
        energy_diff_max = max(diff(log_data$energy_on_rho))*water_density,
        # To simplify the energy interpretation, it's easier to look at the model
        # prior to boundary fluxes occurring. 
        energy_diff_max_while_closed = max(diff(log_data$energy_on_rho[1:before_bc]))*water_density,
        # To further simplify interpretation of energy fluxes, it's nice to
        # measure it relative to the energy. This is because we expect small
        # numerical energy increases [e.g. 1/1000] are no problem -- they seem
        # common with the linear solver in the first few timesteps.
        energy_diff_max_while_closed_relative = max(
            diff(log_data$energy_on_rho[1:before_bc])/log_data$energy_on_rho[1:(before_bc-1)]),
        energy_max_while_closed = max(log_data$energy_on_rho[1:before_bc]*water_density),
        # Mass conservation error               
        mass_error = max(abs(log_data$mass_error)),
        log_file = log_data$log_file
        )

}

# Log files for MSL=0.8 runs, and MSL=0 runs, ordered
all_log_files = c(
    Sys.glob('../../swals/OUTPUTS/ptha18_tonga_MSL0.8/*/RUN*/multidomain*0001.log'),
    Sys.glob('../../swals/OUTPUTS/ptha18_tonga_MSL0/*/RUN*/multidomain*0001.log'))

library(parallel)
MC_CORES = 48
check_data = mclapply(all_log_files, f<-function(logfile){
    tmp=get_time_energy_stagerange_from_log(logfile)        
    output = check_log(tmp)
    rm(tmp)
    return(output)
    }, mc.cores=MC_CORES)

all_checks = do.call(rbind, check_data)

all_checks$ic_name = substring(basename(dirname(dirname(all_log_files))), 1, 59)
# Confirm the runs with different MSL line-up OK
stopifnot(all(all_checks$ic_name[1:1101] == all_checks$ic_name[1102:(1102+1100)]))

# Initially one of the runs with MSL=0.8 died -- all the others finished.
# The problematic run was sorted out as discussed in ../../swals/README.md,
# so we end up with no unfinished runs.
summary(all_checks$run_finished)
#   Mode    TRUE 
#logical    2202 

## First time, we got a failure.
## Here is the failed-run's logfile name
#all_checks$log_file[!all_checks$run_finished]
## [1] "../../swals/OUTPUTS/ptha18_tonga_MSL0.8/ptha18_random_scenarios_kermadectonga2_row_0043831_Mw_95_HS-risetime_0-ambientsealevel_0.8-full-linear_with_manning-0.035-highres_tonga/RUN_20201122_121533834/multidomain_log_image_00000000000000000001.log"
## From inspection this results in a depth raster that is NA everywhere.
## That means it gets treated as having a zero depth in the raster based analysis. The event is very extreme and will have a very small part of the rate, so the effects are small.

# The mass conservation errors are very small and consistent with
# double-precision round-off error (noting the full domain has a volume around
# 2.18e+17),
summary(all_checks$mass_error)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  14.68   24.15   27.66   64.05   32.04  681.27 

# Before the boundary fluxes begin, any energy increases are tiny. It is common to see 'very slightly'
# higher energy just after the start of the simulation, which I think just reflects that the algorithms
# used herein are not precisely energy conservative.
summary((all_checks$energy_max_while_closed - all_checks$energy_start)/all_checks$energy_start)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 0.000e+00 6.516e-06 1.387e-04 1.305e-04 3.132e-03 


# After the boundary fluxes begin it is less straightforward to interpret the potential energy and
# energy, especially when MSL is nonzero. In this plot we compare the energy in runs with MSL=0 and
# MSL=0.8, based on the expected relationships discussed in comments at the start of this script.
# We find that the energies in both models are very similar once accounting for these factors, as expected.
# Note we do not necessarily expect them to be identical [because changing the MSL will change the effect of topography, etc,
# it is not simply an identical model with a different datum, although at global scales we expect them to be very similar].
pdf('Energy_comparison_models_same_source_different_MSL.pdf', width=15, height=7.5)
N = length(all_log_files)/2
for(i in 1:N){
    print(i)
    plot_panel<-function(i){
        # Relation between the energies with MSL=0 and MSL=0.8
        t0 = get_time_energy_stagerange_from_log(all_log_files[i+N])
        t8 = get_time_energy_stagerange_from_log(all_log_files[i])

        # Check the sources are identical
        match_OK = (all_checks$ic_name[i] == all_checks$ic_name[i+N])

        if(any(!is.finite(t8$normalised_energy_on_rho)) |
           any(!is.finite(t0$normalised_energy_on_rho))){
            par(mfrow=c(1,1))
            plot(c(0, 1), c(0, 1), main='Error in this panel due to NA energy')
            # Deliberately throw a try-error
            match_OK = try(log('a'))
        }else{

            # 
            # Compare the energies of the models with MSL=0 and MSL=0.8 
            # These should not necessarily be identical, but they are expected
            # to be pretty similar.
            #
            par(mfrow=c(1,2))
            plot(t8$normalised_energy_on_rho - t8$normalised_energy_on_rho[1], 
                 t0$normalised_energy_on_rho - t0$normalised_energy_on_rho[1],
                 main='Normalised energy in each model \n (accounting for volume change and different MSL)')
            abline(0, 1, col='red')
            plot(t8$kinetic_energy_on_rho[-1], t0$kinetic_energy_on_rho[-1],
                 main=paste0('Kinetic energy in each model \n', all_checks$ic_name[i]))
            abline(0, 1, col='red')
        }
        return(match_OK)
    }

    match_OK = try(plot_panel(i))
    # Report issues
    if(is(match_OK, 'try-error')){
       print(paste0('Try error on i=', i))
    }else if(!match_OK){
       stop(paste0('Failed match on i=', i))
    }
}
dev.off()

