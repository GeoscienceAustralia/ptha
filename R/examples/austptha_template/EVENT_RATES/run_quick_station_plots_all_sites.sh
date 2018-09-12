source R_modules.sh
for i in {1..16}; do
    Rscript quick_station_plots_all_sites.R 0 2 $i 16 > /dev/null &
done
#wait
