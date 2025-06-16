output_dirs <- Sys.glob("ptha_highres_domains_max_stage_at_epistemic_uncertainty_percentile_0.84_exrate_0.0004/")

domains_to_cover <- seq(55, 133, by = 1)

missing_dominains <- list()

for (output_dir in output_dirs) {
    tif_files <- Sys.glob(paste0(output_dir, "/*.tif"))
    domain_numbers <- as.numeric(gsub(".*_domain_index_", "", gsub(".tif", "", tif_files)))
    missing_dominains[[output_dir]] <- setdiff(domains_to_cover, domain_numbers)
}
print(missing_dominains)
