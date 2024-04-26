# Make list of breakwalls. 
# DELIBERATELY IGNORE BUNBURY FLOOD-GATE - the code will manually add it later.
all_breakwalls = c(Sys.glob('perth/*.csv'), Sys.glob('bunbury_busselton/*.csv'), Sys.glob('wonnerup_floodgate/*.csv'))
cat(all_breakwalls, file='swals_breakwall_files.txt', sep="\n")

# This version includes updates for the Vasse Diversion drain (as well as the
# older ones -- that works because the walls are burned into the grid using a "max" operation).
all_breakwalls = c(all_breakwalls, Sys.glob('Vasse_diversion_drain_wall_update/VasseDiversionDrainWallPatch*.csv'),
    Sys.glob('PortGeographUpdate/portGeographe*.csv'))
cat(all_breakwalls, file='swals_breakwall_files_updates2024.txt', sep="\n")

