get_deletion_times_SWALS<-function(multidomain_dir, deletion_policy_days=90){
    x = multidomain_dir
    x_end = substring(x, nchar(x)-17, nchar(x)-3)
    creation_dates = strptime(x_end, format='%Y%m%d_%H%M%S')
    # Deleted after 90 days
    deletion_time = creation_dates - Sys.time() + 
        as.difftime(deletion_policy_days, units = 'days')

    output = data.frame(dir=multidomain_dir, creation_time=creation_dates, 
        time_to_delete=deletion_time)
    return(output)
}


 
