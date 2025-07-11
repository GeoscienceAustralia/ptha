# Code to check the log files, and notes about the interpretation of energy. Run like:
#   Rscript check_log_files.R msl_linear expected_finish_time_hours "../path/that/globs*/log/files/like/multidomain*00001.log"
#

# 
# NOTES ON INTERPRETATION OF SWALS ENERGY CALCULATIONS
#
# As outlined below, the SWALS energy calculations can be affected by the following:
# - Boundary fluxes -- in general energy conservation type results apply to closed domains,
#   whereas in this model the domain is only closed until the tsunami reaches the model boundary,
#   at which point the boundaries are radiating. The 'available potential energy' does TEND to 
#   decrease due to boundary radiation, but there is a non-trivial relationship between this and
#   the SWALS energy (unless we can neglect wetting and drying AND MSL=0).
# - Okada initial conditions with non-zero volume -- this will make the SWALS
#   potential energy differ from the available potential energy when (MSL_LINEAR != 0)
# - Wetting and drying -- SWALS potential energy computations are valid in this
#   case, but the (sum of the kinetic and available potential energy) / rho
#   is no longer a conserved quantity (as commented in
#   Davies et al., 2020, Global Dissipation Models ..., Frontiers in Earth
#   Science).
# - If we run a model at different mean-sea-levels with the same initial
#   conditions, then it is possible to derive the expected relation between the
#   energy in the models (ignoring differences caused by the different MSL). We
#   get an equation like:
#       normalised_energy_on_rho = energy_on_rho + grav * msl_linear * (md_volume[1] - md_volume)
#   where energy_on_rho is what SWALS_reports, md_volume is the time-varying volume in the multidomain, 
#   and msl_linear is the SWALS msl_linear. This is a good way to compare the energy in models that
#   are "the same" except for different MSL_LINEAR values.
# - A simpler comparison is to look at the kinetic energies! 
#
########
#
# THEORY
#
########
#
# Below are notes on the relation between the SWALS potential energy term, and the "available potential energy".
# We generally work with those quantities divided by the density "rho". 
#
# In SWALS we compute potential energy as:
#    SWALS_potential_energy_on_rho = g/2 * integral{ stage^2 - elev^2 } - CONSTANT
# where:
#    CONSTANT = g/2 * integral_where_elevation_less_than_MSL_LINEAR { MSL_LINEAR^2 - elev^2 }
# is designed to make the potential energy zero when the wet-stage is equal to MSL_LINEAR everywhere.
#
# How does this compare with the classical equation for "available potential energy" which
# ignores wetting and drying?:
#    available_potential_energy_on_rho = g/2 * integral_in_wet_areas{ (stage - MSL_LINEAR)^2 }
# 
# Supposing there is no wetting and drying, we can expand the SWALS potential energy as:
#     SWALS_potential_energy_on_rho = g/2 * integral_in_wet_areas{ stage^2 - MSL_LINEAR^2}
# whereas the classical equation expands as.
#     available_potential_energy_on_rho = g/2 * integral_in_wet_areas{ stage^2 - 2*stage * MSL_LINEAR + MSL_LINEAR^2}
# If we suppose the mean_stage REALLY is MSL_LINEAR, i.e.
#    integral_in_wet_areas{stage} == integral_in_wet_areas{MSL_LINEAR}
# then 
#    integral_in_wet_areas{ -2*stage*MSL_LINEAR } = integral_in_wet_areas{ -2 * MSL_LINEAR^2 }
# so we can rearrange the classical equation as 
#     available_potential_energy_on_rho = g/2 * integral_in_wet_areas{ stage^2 - MSL_LINEAR^2} = SWALS_potential_energy_on_rho
# In summary, the SWALS calculation is the same as the classical available potential energy if there
# is no wetting and drying, and the mean-wet-stage == MSL_LINEAR.
#
# In general that might not be the case, however [e.g. if we initialise the
# model with an Okada initial condition that has non-zero volume, the
# mean-wet-stage will no longer be equal to the models MSL]. In that case, if
# there is no wetting and drying, then [KEY EQUATION]:
#     available_potential_energy_on_rho = 
#         SWALS_potential_energy_on_rho + g/2 * 2 * MSL_LINEAR*integral_in_wet_areas{ (MSL_LINEAR - stage) } = 
#         SWALS_potential_energy_on_rho - g*MSL_LINEAR*{ CHANGE_IN_VOLUME + SOME_CONSTANT }
# Notice that the time-variation of the adjustment term is related to the
# volume of water in the domain, which is something our volume conservation
# analysis keeps track of. Also, the latter term vanishes if MSL_LINEAR=0,
# which greatly simplifies that case.
#
# So with MSL not equal to 0, the SWALS_potential_energy_on_rho may be greater
# than or less than the available_potential_energy_on_rho. We should still have
# conservation of the SWALS energy in closed domains [ as the SWALS
# potential energy is the actual potential energy plus some constant], or
# dissipation in models with friction, but the SWALS potential energy does not
# correspond to the available potential energy (unless MSL_LINEAR=0).  However,
# if we understand the change in volume of the domain over time, and know
# MSL_LINEAR, then we can account for those factors [as done in the function
# below]
#
# Summary: Wetting and drying will invalidate the classical
# available_potential_energy calculation, but not SWALS calculation. However,
# if the mean stage is not MSL_LINEAR, then the SWALS calculation will not
# correspond to the available potential energy. Boundary fluxes make all
# reasoning about the energy balance more tricky. Thus a good approach is to:
# 1) Confirm the energy is non-increasing while the domain is closed.
# 2) When the domain starts radiating, we expect energy decreases for MSL_LINEAR=0,
#    but perhaps not for other values of MSL_LINEAR because the energy is affected
#    by the volume of fluid in the domain. 
# 3) However, we can cross-check the cases with MSL_LINEAR != 0, by comparing with
#    the cases with MSL_linear=0 (which we already checked are well behaved), using
#    the relations above.
# 
library(stringr)

get_time_energy_stagerange_from_log<-function(log_file, msl_linear = 0.0, flux_threshold = 1.0){

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
    
    max_speed_data = as.numeric(x[k-6])

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
    #normalised_energy_on_rho = energy_on_rho + grav * msl_linear * (md_volume[1] - md_volume)
    normalised_energy_on_rho = energy_on_rho + grav * msl_linear * (volume_change[1] - volume_change) # Less roundoff

    
    # Ad-hoc indicator for the time when boundary-fluxes are zero [i.e.  before
    # the wave reaches a model boundary]. After the model starts having
    # boundary fluxes, we can't guarentee anything about energy conservation.
    before_bc = min(which(abs(boundary_flux_integral) > flux_threshold)) - 1
    if(!is.finite(before_bc)) before_bc = length(time)
    before_bc_limit = before_bc
    if(before_bc <= 1) before_bc = 3


    return(list(time=time, 
                energy_on_rho=energy_on_rho, 
                kinetic_energy_on_rho=kinetic_energy_on_rho,
                normalised_energy_on_rho = normalised_energy_on_rho,
                max_stage = max_stage_data, 
                min_stage=min_stage_data,
                max_speed = max_speed_data, 
                log_file=log_file, 
                mass_error = mass_error, 
                boundary_flux_integral = boundary_flux_integral,
                md_volume = md_volume, 
                msl_linear=msl_linear,
                before_bc = before_bc,
                before_bc_limit = before_bc_limit,
                before_bc_time = time[before_bc]))
}

# Logical checks on the result of get_time_energy_stagerange_from_log
check_log<-function(log_data, expected_final_time){

    n = length(log_data$time)

    water_density = 1024

    # Indices related to the time at which boundary radiation seems to begin (based on a ad-hoc flux_threshold)
    before_bc = log_data$before_bc
    before_bc_limit = log_data$before_bc_limit

    data.frame(
        n = n,
        msl_linear = log_data$msl_linear, 
        before_bc = before_bc_limit, # Number of time-steps before the stage-perturbation reaches the boundary, with min=2
        run_finished = (abs(log_data$time[n] - (expected_final_time * 3600)) < 10), # Ran for N hours (last output within 10 seconds of that)
        before_bc_time = log_data$time[before_bc], # Roughly the time at which boundary fluxes begin.
        stage_maxima = max(log_data$max_stage),
        stage_minima = min(log_data$min_stage),
        initial_volume = log_data$md_volume[1],
        energy_start = log_data$energy_on_rho[1] * water_density,
        energy_end = log_data$energy_on_rho[n] * water_density,
        energy_max = max(log_data$energy_on_rho) * water_density,
        energy_max_index = which.max(log_data$energy_on_rho),
        kinetic_energy_max = max(log_data$kinetic_energy_on_rho) * water_density, 
        kinetic_energy_max_index = which.max(log_data$kinetic_energy_on_rho),
        energy_max_while_closed = max(log_data$energy_on_rho[1:before_bc]) * water_density,
        # The change in energy is harder to interpret once boundary fluxes can happen
        energy_diff_max = max(diff(log_data$energy_on_rho))*water_density,
        energy_diff_max_relative = max(
            diff(log_data$energy_on_rho[1:n])/log_data$energy_on_rho[1:(n-1)]),
        # To simplify the energy interpretation, it's easier to look at the model
        # prior to boundary fluxes occurring. 
        energy_diff_max_while_closed = max(diff(log_data$energy_on_rho[1:before_bc]))*water_density,
        # Have a look at what happens once the boundary conditions are radiating
        normalised_energy_diff_max = max(diff(log_data$normalised_energy_on_rho))*water_density,
        # Mass conservation error               
        mass_error = max(abs(log_data$mass_error)),
        max_abs_boundary_flux_integral = max(abs(log_data$boundary_flux_integral)),
        log_file = log_data$log_file
        )

}


MSL_LINEAR_DEFAULT = as.numeric(commandArgs(trailingOnly=TRUE)[1])
EXPECTED_FINAL_TIME = as.numeric(commandArgs(trailingOnly=TRUE)[2])
all_log_files = unique(
    sapply(commandArgs(trailingOnly=TRUE)[-c(1,2)], Sys.glob, USE.NAMES=FALSE))

stopifnot(is.finite(MSL_LINEAR_DEFAULT) & is.finite(EXPECTED_FINAL_TIME))

library(parallel)
MC_CORES = 1 #48
check_data = mclapply(all_log_files, function(logfile){
    tmp=get_time_energy_stagerange_from_log(logfile, msl_linear = MSL_LINEAR_DEFAULT)
    output = check_log(tmp, expected_final_time=EXPECTED_FINAL_TIME)
    rm(tmp)
    return(output)
    }, mc.cores=MC_CORES)

all_checks = do.call(rbind, check_data)

all_checks$ic_name = substring(basename(dirname(all_log_files)), 1, 59)

# Confirm that all the model runs finished
print('Did the model runs finish?')
print(summary(all_checks$run_finished))
## Shut bunbury floodgate
#   Mode    TRUE 
# logical     369 

## Open bunbury floodgate
#   Mode    TRUE 
# logical     369 


# Check the mass conservation errors are very small and consistent with
# double-precision round-off error
print('Mass conservation errors relative to initial volume')
print(summary(all_checks$mass_error/all_checks$initial_volume))
## Shut bunbury floodgate
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 7.222e-17 8.072e-17 8.431e-17 8.794e-17 9.081e-17 1.646e-16 

## Open bunbury Floodgate
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
# 6.899e-17 8.073e-17 8.532e-17 8.838e-17 9.231e-17 1.640e-16 



# Here we compare the maximum mass error with the maximum abs(boundary_flux_integral),
# confirming it is always a negligible fraction of the changes in mass due to fluxes through
# the boundaries.
print('Mass conservation errors relative to boundary flux integral')
print(summary(all_checks$mass_error/all_checks$max_abs_boundary_flux_integral))
## Shut bunbury floodgate
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#3.680e-11 7.870e-11 1.578e-10 8.975e-09 7.981e-10 4.958e-07 

## Open bunbury floodgate
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#3.750e-11 7.690e-11 1.703e-10 9.025e-09 7.780e-10 4.978e-07


# Theoretically we expect energy to decay over time when there are no fluxes in
# or out of the domain. Numerically, below we confirm that before the boundary fluxes 
# begin any energy increases are tiny. It is common to see 'very slightly'
# higher energy just after the start of the simulation, which just reflects
# that the algorithms used herein are not precisely energy conservative. But that
# should be small relative to the available energy (estimated here as 2x'max-kinetic-energy')
print('(Maximum energy - initial energy) relative to (2x maximum kinetic energy), BEFORE BOUNDARY FLUXES')
print('(typically very small unless the source is time-varying)')
print(summary((all_checks$energy_max_while_closed - all_checks$energy_start)/(2*all_checks$kinetic_energy_max)))
## Shut bunbury floodgate
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.000e+00 5.272e-05 2.966e-04 2.843e-04 4.598e-04 7.459e-04

## Open bunbury floodgate
#    Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#0.000e+00 5.272e-05 2.966e-04 2.843e-04 4.598e-04 7.459e-04


# We expect the maximum kinetic energy to occur near the start of the simulation (usually second time-step)
print('Time index with largest kinetic energy (usually near start)')
print(summary(all_checks$kinetic_energy_max_index))
## Shut bunbury floodgate
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  2.000   2.000   3.000   3.477   4.000  11.000

## Open bunbury floodgate
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  2.000   2.000   3.000   3.477   4.000  11.000


# When there are fluxes in and out of the domain, we lose the typical energy
# conservation relation.  But given we have a radiation boundary condition, it
# seems reasonable that the available potential energy should decrease when the
# boundary is radiating. When MSL=0, the SWALS potential energy is the same as the
# available potential energy if wetting and drying is unimportant (this is not 
# true for runs with nonzero MSL, as it also depends on the volume change over
# time). 
print('Maximum energy increase between timesteps, relative to (2x maximum kinetic energy)')
print(', normalised for volume change but not wetting and drying.')
print('(typically very small unless the source is time-varying)')
print(summary(all_checks$normalised_energy_diff_max/(2*all_checks$kinetic_energy_max)))
## Shut bunbury floodgate
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -0.0003201  0.0001092  0.0003299  0.0003299  0.0004850  0.0035417

## Open bunbury floodgate
#      Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
#-0.0003201  0.0001092  0.0003299  0.0003297  0.0004850  0.0038261




# Plot
pdf('Energy_time_series_new.pdf', width=15, height=7.5)
N = length(all_log_files)
for(i in 1:N){
    print(i)
    plot_panel<-function(i, MSL=MSL_LINEAR_DEFAULT){

        t8 = get_time_energy_stagerange_from_log(all_log_files[i], msl_linear = MSL)
        #t0 = get_time_energy_stagerange_from_log(
        #    gsub("sealevel60cm-revised", 'sealevel60cm', all_log_files[i]), msl_linear=0.6)
        #t0 = get_time_energy_stagerange_from_log(
        #    #gsub("sealevel60cm-lowres", 'sealevel60cm', all_log_files[i]), msl_linear=0.6)
        #    gsub("sealevel60cm-reviseddomain-highres", 'sealevel60cm', all_log_files[i]), msl_linear=0.6)

        if(any(!is.finite(t8$normalised_energy_on_rho))){
            par(mfrow=c(1,1))
            plot(c(0, 1), c(0, 1), main='Error in this panel due to NA energy')
            # Deliberately throw a try-error
            match_OK = try(log('a'))
        }else{
            match_OK = TRUE
            # 
            # Compare the energies of the models with MSL=0 and MSL=0.8 
            # These should not necessarily be identical, but they are expected
            # to be pretty similar.
            #
            model_file = all_log_files[i] # all_checks$ic_name[i]
            model_id = str_extract(model_file, "(?<=_row_).*?(?=_Mw_)")
            model_mw = as.double(str_extract(model_file, "(?<=_Mw_).*?(?=_HS-)"))/10
            model_name = paste0("PTHA18 ID: ", model_id, " Mw=", model_mw)

            par(mfrow=c(2,2))
            time_h = t8$time / 3600
            t8$before_bc_time = t8$before_bc_time / 3600
            plot(
                time_h,
                t8$normalised_energy_on_rho,
                t='l',
                main=paste0('Normalised Energy (accounting for volume change) \n', model_name),
                xlab='Time (h)',
                ylab="Normalised Energy (J m^3 / kg)"
            )
            abline(v=t8$before_bc_time, col='blue', lty='dashed')
            #points(t0$time, t0$normalised_energy_on_rho, t='l', col='blue')
            points(
                time_h,
                t8$energy_on_rho,
                t='l',
                col='red',
                lty='dashed'
            )
            legend(
                'topright',
                c('Energy', 'Energy normalised for volume change', 'Time when |boundary fluxes| > threshold'), 
                col=c('red', 'black', 'blue'),
                lty=c('dashed', 'solid', 'dashed'),
                bty='n'
            )

            plot(
                time_h,
                t8$kinetic_energy_on_rho,
                t='l', 
                main=paste0('Kinetic energy in model \n', model_name),
                xlab='Time (h)',
                ylab="Normalised Kinetic Energy (J m^3 / kg)"
            )
            abline(v=t8$before_bc_time, col='blue', lty='dashed')
            #points(t0$time, t0$kinetic_energy_on_rho, t='l', col='blue')
            plot(
                time_h,
                t8$max_stage,
                t='l',
                ylim=c(0, max(t8$max_stage)), 
                main=paste0('Max-stage in model \n', model_name),
                xlab='Time (h)',
                ylab='Max-stage (m)'
            )
            abline(v=t8$before_bc_time, col='blue', lty='dashed')
            #points(t0$time, t0$max_stage, t='l', col='blue')

            plot(
                time_h,
                t8$max_speed,
                t='l',
                ylim=c(0, max(t8$max_speed)), 
                main=paste0('Max-speed in model \n', model_name),
                xlab='Time (h)',
                ylab='Max-speed (m/s)'
            )
            abline(v=t8$before_bc_time, col='blue', lty='dashed')
            #points(t0$time, t0$max_speed, t='l', col='blue')
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
