#
# Plot the model-vs-data at Nukualofa for validation runs.
#
 
source('plots/plot_all.R', chdir=TRUE)
all_output_dirs = Sys.glob('./OUTPUTS/*/RUN*')

output_dir = 'plots/historic_events_time_series_plots/'
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)

# Find the directories with the validation runs
all_validation_info = lapply(
    c('Tonga2006', 'Chile2010', 'Tohoku2011', 'Chile2015'),   
    f<-function(x){
       dirnames = all_output_dirs[grep(paste0(x, '_validation_'), all_output_dirs, fixed=TRUE)]
       event_name = rep(x, length=length(dirnames))
       return(list(dirnames=dirnames, event_name=event_name))
    })
all_validation_md_dir     = unlist(lapply(all_validation_info, f<-function(x) x$dirnames))
all_validation_event_name = unlist(lapply(all_validation_info, f<-function(x) x$event_name))

# For each directory, make the plot
for(i in 1:length(all_validation_md_dir)){
    md_dir = all_validation_md_dir[i]
    event_name = all_validation_event_name[i]

    plot_output_dir = paste0(output_dir, event_name)
    dir.create(plot_output_dir, showWarnings=FALSE)

    plot_file_name = paste0(plot_output_dir, '/', 'nukualofa_gauge_modelVdata_',
        gsub('/', '-', paste0(basename(dirname(md_dir)), '-', basename(md_dir))), 
        '.png')
    # Make a plot inside the multidomain_directory of each validation run 
    png(plot_file_name, width=10, height=4, units='in', res=300)
    historic_event_gauge_plot(md_dir)
    dev.off()
}

