output_dirs <- Sys.glob("ptha/sea_level_vary/highres_depth_with_variance/*")

domains_to_cover <- seq(55, 133, by = 1)

missing_dominains <- list()

for (output_dir in output_dirs) {
    tif_files <- Sys.glob(paste0(output_dir, "/*.tif"))
    domain_numbers <- as.numeric(gsub(".logic_tree_mean_HS_domain_2_depth_as_max_stage_minus_elevation0_domain__exceedance_rate_with_threshold_0.001", "", gsub(".tif", "", tif_files)))
    missing_dominains[[output_dir]] <- setdiff(domains_to_cover, domain_numbers)
}
print(missing_dominains)
