# Choose sites that (roughly) match the different static sea levels
get_sites_and_sea_levels<-function(){

    sites_at_sea_levels = list(
        'Geraldton'  = list(sealevel=0.65, plot_bbox=c(114.5633 , -28.7998 , 114.6169 , -28.7611),  domain_indices = c(267, 269)),
        'Kalbarri'   = list(sealevel=0.65, plot_bbox=c(114.14155, -27.72868, 114.19480, -27.69112), domain_indices=274),
        'ShelterBay' = list(sealevel=0.65, plot_bbox=c(113.18655, -26.185, 113.21325, -26.1587),  domain_indices=c(279)),
        'Carnarvon'  = list(sealevel=0.95, plot_bbox=c(113.6170 , -24.910  , 113.6784 , -24.8663),  domain_indices=seq(282, 285)),
        'CoralBay'   = list(sealevel=0.95, plot_bbox=c(113.7603 , -23.1506 , 113.7834 , -23.1341),  domain_indices=c(320, 321)),
        'Exmouth'    = list(sealevel=1.17, plot_bbox=c(114.112  , -21.9615 , 114.1593 , -21.92788), domain_indices=c(286)),
        'Onslow'     = list(sealevel=1.42, plot_bbox=c(115.094  , -21.66   , 115.1407 , -21.6263),  domain_indices = c(290, 291))
    )

    # Check bbox is ordered OK
    stopifnot(all(unlist(lapply(sites_at_sea_levels, 
        function(x) (x$plot_bbox[1]<x$plot_bbox[3]) & (x$plot_bbox[2]<x$plot_bbox[4]) ))))

    sites_at_sea_levels
}
sites_at_sea_levels = get_sites_and_sea_levels()


