library(ggplot2)


plot_h95 <- function(model_file, ATWS_ZONE_NAME, T2_file) {
    model_H = read.csv(model_file)
    t2_H = read.csv(T2_file)
    merged_data <- merge(model_H, t2_H, by='event_name', suffixes=c("_model","_t2"))
    print(merged_data)
    print(model_H)
    print(t2_H)

    # Plot the JATWC_H values to compare model against t2 values for each corresponding event
    # like a predicted (model) vs observed (t2) plot
    p <- ggplot(merged_data, aes(x=jatwc_H_model, y=jatwc_H_t2, color=event_name)) +
        geom_point(size=2) +
        geom_abline(slope=1, intercept=0, linetype="dashed", color="black") +
        labs(
        x = "Model H95",
        y = "T2 H95",
        color = "Event",
        title = paste("Model vs T2", ATWS_ZONE_NAME)
        )
    
    ggsave(paste0('JATWC_H95_comparison_', ATWS_ZONE_NAME, '.png'), plot=p, width=8, height=6, units='in')
}

model_file = 'all-JATWC-H_Capricornia-Coast.csv'
t2_file = 'JATWC-H95-T2.csv'
ATWS_ZONE_NAME = 'Capricornia Coast'
plot_h95(model_file, ATWS_ZONE_NAME, t2_file)