# Make list of breakwalls. 
# DELIBERATELY IGNORE BUNBURY FLOOD-GATE - the code will manually add it later.
all_breakwalls = c(Sys.glob('perth/*.csv'), Sys.glob('bunbury_busselton/*.csv'), Sys.glob('wonnerup_floodgate/*.csv'))
cat(all_breakwalls, file='swals_breakwall_files.txt', sep="\n")
