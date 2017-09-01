all_Rdata = Sys.glob('../SOURCE_ZONES/*/TSUNAMI_EVENTS/plots/*.Rdata')

uniform_store = list()
stochastic_store = list()
variable_uniform_store = list()

for(i in 1:length(all_Rdata)){

    print(all_Rdata[i])
    event_env = new.env()

    load(all_Rdata[i], envir=event_env)

    ngauges = length(event_env$uniform_slip_stats)

    uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )
    stochastic_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$stochastic_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )
    variable_uniform_store[[i]] = list(
        data= matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        model=matrix(NA, ncol=ngauges, nrow=length(event_env$variable_uniform_slip_stats[[1]])),
        nunique_models = rep(NA, ngauges)
        )

    for(gauge_ind in 1:ngauges){

        # Uniform
        uniform_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$uniform_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        uniform_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$uniform_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))

        uniform_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(uniform_store[[i]]$model[,gauge_ind]))

        # Stochastic
        stochastic_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$stochastic_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        stochastic_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$stochastic_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))
        stochastic_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(stochastic_store[[i]]$model[,gauge_ind]))

        # variable_uniform
        variable_uniform_store[[i]]$data[,gauge_ind] = 
            unlist(lapply(event_env$variable_uniform_slip_stats[[gauge_ind]],
                f<-function(x) diff(x$data_range)))
        variable_uniform_store[[i]]$model[,gauge_ind] = 
            unlist(lapply(event_env$variable_uniform_slip_stats[[gauge_ind]], 
                f<-function(x) diff(x$model_range)))
        variable_uniform_store[[i]]$nunique_models[gauge_ind] = 
            length(unique(variable_uniform_store[[i]]$model[,gauge_ind]))
        
    }

}

#
# Compute the fraction of model events that had stage-range exceeding the data.
# If multiple gauges exist, take the median value over all gauges.
#
meds_uniform          = unlist(lapply(uniform_store         , f<-function(x) median(colMeans(x$model > x$data))))
# [1] 0.4166667 0.8333333 0.0000000 0.0000000 0.9285714 0.1785714 0.3157895 0.8333333 0.3484848 0.0000000 0.5000000 0.5714286
meds_stochastic       = unlist(lapply(stochastic_store      , f<-function(x) median(colMeans(x$model > x$data))))
# [1] 0.6444444 0.8944444 0.6176471 0.1111111 0.7952381 0.1761905 0.4438596 0.8244444 0.4939394 0.2825397 0.7595238 0.5833333
meds_variable_uniform = unlist(lapply(variable_uniform_store, f<-function(x) median(colMeans(x$model > x$data))))
# [1] 0.35555556 0.65555556 0.18823529 0.02222222 0.49047619 0.11904762 0.25964912 0.61666667 0.27474747 0.14920635 0.50714286 0.44047619

meds2_uniform          = (lapply(uniform_store         , f<-function(x) (colMeans(x$model > x$data))))
meds2_stochastic       = (lapply(stochastic_store      , f<-function(x) (colMeans(x$model > x$data))))
meds2_variable_uniform = (lapply(variable_uniform_store, f<-function(x) (colMeans(x$model > x$data))))

## > meds2_uniform
## [[1]]
## [1] 0.7500000 0.6666667 0.4166667 0.2500000 0.0000000
## 
## [[2]]
## [1] 0.8333333
## 
## [[3]]
##  [1] 0.11764706 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.23529412 0.00000000 0.05882353 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000 0.00000000
## [17] 0.00000000 0.00000000 0.00000000 0.00000000 0.58823529 0.29411765 0.00000000 0.29411765 0.35294118 0.41176471 0.52941176 0.64705882
## 
## [[4]]
## [1] 0 0 0 0 0
## 
## [[5]]
## [1] 0.9285714
## 
## [[6]]
## [1] 0.0000000 0.3571429
## 
## [[7]]
## [1] 0.05263158 0.50000000 0.13157895 0.00000000 0.34210526 0.44736842 0.21052632 0.31578947 0.31578947
## 
## [[8]]
## [1] 0.7333333 0.8666667 0.8000000 0.9666667
## 
## [[9]]
##  [1] 0.93939394 0.27272727 0.12121212 0.24242424 0.33333333 0.06060606 0.72727273 0.33333333 0.66666667 0.24242424 0.48484848 0.36363636 0.36363636 0.51515152 0.39393939 0.30303030
## 
## [[10]]
## [1] 0
## 
## [[11]]
## [1] 0.5
## 
## [[12]]
## [1] 0.8571429 0.2857143
## 
## > meds2_stochastic
## [[1]]
## [1] 0.8444444 0.7944444 0.6444444 0.4944444 0.2166667
## 
## [[2]]
## [1] 0.8944444
## 
## [[3]]
##  [1] 0.6823529 0.6588235 0.5333333 0.3921569 0.3686275 0.5490196 0.8980392 0.6235294 0.8078431 0.5058824 0.6705882 0.3529412 0.6039216 0.6117647 0.5137255 0.5568627 0.6862745
## [18] 0.6000000 0.2784314 0.4392157 0.9254902 0.7607843 0.4901961 0.8352941 0.8901961 0.8470588 0.9607843 0.9333333
## 
## [[4]]
## [1] 0.22222222 0.06666667 0.13333333 0.11111111 0.04444444
## 
## [[5]]
## [1] 0.7952381
## 
## [[6]]
## [1] 0.05714286 0.29523810
## 
## [[7]]
## [1] 0.4280702 0.6789474 0.4245614 0.2789474 0.4982456 0.4438596 0.4140351 0.6421053 0.5421053
## 
## [[8]]
## [1] 0.8044444 0.9133333 0.8444444 0.8022222
## 
## [[9]]
##  [1] 0.8606061 0.5595960 0.4545455 0.5838384 0.6929293 0.5313131 0.7070707 0.4949495 0.6828283 0.3191919 0.4323232 0.4444444 0.4181818 0.4929293 0.4646465 0.4484848
## 
## [[10]]
## [1] 0.2825397
## 
## [[11]]
## [1] 0.7595238
## 
## [[12]]
## [1] 0.7952381 0.3714286
## 
## > meds2_variable_uniform
## [[1]]
## [1] 0.71111111 0.58888889 0.35555556 0.20555556 0.06666667
## 
## [[2]]
## [1] 0.6555556
## 
## [[3]]
##  [1] 0.35686275 0.25882353 0.12941176 0.05490196 0.15686275 0.16078431 0.62745098 0.18431373 0.55686275 0.16078431 0.23921569 0.05882353 0.15294118 0.19215686 0.14901961 0.12156863
## [17] 0.24705882 0.18039216 0.04705882 0.11764706 0.64313725 0.50980392 0.12941176 0.41176471 0.42745098 0.40784314 0.75294118 0.78823529
## 
## [[4]]
## [1] 0.04444444 0.02222222 0.04444444 0.02222222 0.00000000
## 
## [[5]]
## [1] 0.4904762
## 
## [[6]]
## [1] 0.04285714 0.19523810
## 
## [[7]]
## [1] 0.1508772 0.5245614 0.1929825 0.1228070 0.3192982 0.2596491 0.2087719 0.3157895 0.3526316
## 
## [[8]]
## [1] 0.5666667 0.8244444 0.5533333 0.6666667
## 
## [[9]]
##  [1] 0.7636364 0.2626263 0.1959596 0.3030303 0.3737374 0.2161616 0.4545455 0.2404040 0.4848485 0.1838384 0.2808081 0.3212121 0.2424242 0.3212121 0.2686869 0.2626263
## 
## [[10]]
## [1] 0.1492063
## 
## [[11]]
## [1] 0.5071429
## 
## [[12]]
## [1] 0.6428571 0.2380952

save.image('model_data_envelope_summary_statistics.Rdata')

#
# Info on number of unique models [since double-ups can occur, e.g.
# if 2 variable_uniform slip models have the same area and Mw. This
# even happens for single unit-source stochastic models
#
for(i in 1:12){
    print(paste0('################', i))
    print(basename(all_Rdata[i]))
    print(c(uniform_store[[i]]$nunique_models, '/', length(uniform_store[[i]]$model[,1])))
    print(c(stochastic_store[[i]]$nunique_models, '/', length(stochastic_store[[i]]$model[,1])))
    print(c(variable_uniform_store[[i]]$nunique_models, '/', length(variable_uniform_store[[i]]$model[,1])))
}

## [1] "################1"
## [1] "gauge_summary_stats_session_kermadectonga_samoa_2009_09_29_Mw8.1.Rdata"
## [1] "12" "12" "12" "12" "12" "/"  "12"
## [1] "180" "180" "180" "180" "180" "/"   "180"
## [1] "127" "127" "127" "127" "127" "/"   "180"
## [1] "################2"
## [1] "gauge_summary_stats_session_kermadectonga_tonga_2009_03_19_Mw7.7.Rdata"
## [1] "12" "/"  "12"
## [1] "178" "/"   "180"
## [1] "136" "/"   "180"
## [1] "################3"
## [1] "gauge_summary_stats_session_kurilsjapan_tohoku_2011_03_11_Mw9.1.Rdata"
##  [1] "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17"
## [16] "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "17" "/"  "17"
##  [1] "255" "255" "255" "255" "255" "255" "255" "255" "255" "255" "255" "255"
## [13] "255" "255" "255" "255" "255" "255" "255" "255" "255" "255" "255" "255"
## [25] "255" "255" "255" "255" "/"   "255"
##  [1] "230" "230" "230" "230" "230" "230" "230" "230" "230" "230" "230" "230"
## [13] "230" "230" "230" "230" "230" "230" "230" "230" "230" "230" "230" "230"
## [25] "230" "230" "230" "230" "/"   "255"
## [1] "################4"
## [1] "gauge_summary_stats_session_newhebrides_northNewHebrides_2013_02_06_Mw7.9.Rdata"
## [1] "3" "3" "3" "3" "3" "/" "3"
## [1] "45" "45" "45" "45" "45" "/"  "45"
## [1] "31" "31" "31" "31" "31" "/"  "45"
## [1] "################5"
## [1] "gauge_summary_stats_session_newhebrides_vanuatu_north_2009_10_07_Mw7.8.Rdata"
## [1] "14" "/"  "14"
## [1] "210" "/"   "210"
## [1] "102" "/"   "210"
## [1] "################6"
## [1] "gauge_summary_stats_session_puysegur_puysegur_2009_07_15_Mw78.Rdata"
## [1] "14" "14" "/"  "14"
## [1] "205" "205" "/"   "210"
## [1] "98"  "98"  "/"   "210"
## [1] "################7"
## [1] "gauge_summary_stats_session_southamerica_chile_2010_02_27_Mw8.8.Rdata"
##  [1] "38" "38" "38" "38" "38" "38" "38" "38" "38" "/"  "38"
##  [1] "570" "570" "570" "570" "570" "570" "570" "570" "570" "/"   "570"
##  [1] "509" "509" "509" "509" "509" "509" "509" "509" "509" "/"   "570"
## [1] "################8"
## [1] "gauge_summary_stats_session_southamerica_southamerica_2014_04_01_Mw8.2.Rdata"
## [1] "30" "30" "30" "30" "/"  "30"
## [1] "450" "450" "450" "450" "/"   "450"
## [1] "361" "361" "361" "361" "/"   "450"
## [1] "################9"
## [1] "gauge_summary_stats_session_southamerica_southamerica_2015_09_16_Mw8.3.Rdata"
##  [1] "33" "33" "33" "33" "33" "33" "33" "33" "33" "33" "33" "33" "33" "33" "33"
## [16] "33" "/"  "33"
##  [1] "495" "495" "495" "495" "495" "495" "495" "495" "495" "495" "495" "495"
## [13] "495" "495" "495" "495" "/"   "495"
##  [1] "412" "412" "412" "412" "412" "412" "412" "412" "412" "412" "412" "412"
## [13] "412" "412" "412" "412" "/"   "495"
## [1] "################10"
## [1] "gauge_summary_stats_session_sunda_mentawai_2010_10_25_Mw7.9.Rdata"
## [1] "21" "/"  "21"
## [1] "312" "/"   "315"
## [1] "211" "/"   "315"
## [1] "################11"
## [1] "gauge_summary_stats_session_sunda_sumatra_2007_09_12_Mw8.5.Rdata"
## [1] "28" "/"  "28"
## [1] "420" "/"   "420"
## [1] "391" "/"   "420"
## [1] "################12"
## [1] "gauge_summary_stats_session_sunda_sumatra_2010_04_06_Mw7.8.Rdata"
## [1] "14" "14" "/"  "14"
## [1] "206" "206" "/"   "210"
## [1] "140" "140" "/"   "210"
