#
# Extract spatial coordinates of "good nearshore" gauges in the highres domains
# that are close enough to our points of interest to be included. 
#

event_stats = readRDS('event_stats.RDS')

# Find sites that are actually used in our statistics
k = which(event_stats$good_nearshore & event_stats$is_gauge_in_highres_domain & event_stats$distance_to_gauge < 200)
event_stats_k = event_stats[k,]

#
# Get each "site name" only once. This does not completely avoid double-ups
# because often we use different names for the same site at different times.
#
unique_gauges = unique(event_stats_k$sites)
N = match(unique_gauges, event_stats_k$sites)
output_df = data.frame(sites = event_stats_k$sites[N], lon=event_stats_k$obs_lon[N], lat=event_stats_k$obs_lat[N])
write.csv(output_df, file='gauge_location_info.csv', row.names=FALSE)

#
# Count gauges observations for each event
#
unique_sites_events = unique(event_stats_k$site_and_event)
events_only = unlist(lapply(unique_sites_events, function(x) strsplit(x, split="_")[[1]][1]))
nc1 = nchar(events_only)+2
nc2 = nchar(unique_sites_events)
sites_only = substring(unique_sites_events, nc1, nc2)

print(unique_sites_events)
print(table(events_only))

#
# Clean up the table for easier communication in the paper. We want to show
# - lon, lat, name of a particular gauge site, without the complexities caused
#   by repeated names for the same site (like appears in the underlying database,
#   for good reason).
# - the tsunami events that we use at that site
# - the data provider name
#

# Use to associated a data provider to the data
match_data_source_to_name<-function(site_name){
    # Here site_name is an entry of names(GAUGE_DATA)
    if(grepl('WADOT', site_name) | grepl('WADoT', site_name)){
        return('WADOT')
    }else if(grepl('DPIE', site_name)){
        return('NSW DPIE')
    }else if(grepl('BOM', site_name)){
        return('BOM')
    }else if(grepl('_PA', site_name)){
        return('PANSW')
    }else if(grepl('FremantleInnerHarbour_5min', site_name)){
        return("FP")
    }else if(grepl("GoldCoastSandBypass", site_name)){
        return('BOM')
    }else if(grepl("Sydney_FortDenison1960", site_name)){
        return('W2018')
    }else if(grepl("Cronulla_CSIRO_Fisheries_1960", site_name)){
        return('D2020')
    }else{
        return(NA)
    }
}

# Use to associate a nice name to the data
match_site_to_station_name = list(
    "BallinaBreakwall_1min_DPIE" = 'Ballina Breakwall',
    "Bermagui_1min_DPIE" = 'Bermagui',                    
    "BrunswickHeads_1min_DPIE" = 'Brunswick Heads',               
    "CoffsHarbourInnerPumpoutJetty_1min_DPIE" = 'Coffs Harbour Inner Pumpout Jetty',
    "CrookhavenHeads_1min_DPIE" = 'Crookhaven Heads',              
    "CrowdyHeadFishermansWharf_1min_DPIE" = 'Crowdy Head Fishermans Wharf',    
    "Forster_1min_DPIE" = 'Forster',                      
    "LordHoweIsland_1min_DPIE" = 'Lord Howe Island',               
    "PortMacquarie_1min_DPIE" = 'Port Macquarie',                 
    "ShoalBay_1min_DPIE" = 'Shoal Bay',                     
    "TweedEntranceSouth_1min_DPIE" = 'Tweed Entrance South',           
    "Yamba_1min_DPIE" = 'Yamba',                        
    "Eden_1min_DPIE" = 'Eden',                          
    "TwofoldBay_1min_BOM" = 'Twofold Bay',                    
    "GoldCoastSandBypass" = 'Gold Coast Sand Bypass',                    
    "BatemansBay_PrincessJetty_1min_DPIE" = 'Batemans Bay Princess Jetty',    
    "Ulladulla_1min_DPIE" = 'Ulladulla',                    
    "JervisBay_1min_DPIE" = 'Jervis Bay',                    
    "Bundeena_1min_DPIE" = 'Bundeena',                     
    "Sydney_MiddleHarbour_1min_DPIE" = 'Sydney Middle Harbour',         
    "Sydney_FortDenison1960" = 'Sydney Fort Denison',                 
    "Cronulla_CSIRO_Fisheries_1960" = 'Cronulla',          
    "Sydney_FortDenison_1min_PA" = 'Sydney Fort Denison',             
    "Sydney_FortDenison_1min_PA_b" = 'Sydney Fort Denison',           
    "Sydney_BotanyBay_1min_PA" = 'Botany Bay Bulk Liquids Berth',                
    "Sydney_Botany_Bay_Pilot_Jetty_1min_PA" = 'Botany Bay Pilot Jetty', 
    "Sydney_Botany_Bay_Kurnell_Wharf_1min_PA" = 'Botany Bay Kurnell Wharf',
    "Hawkesbury_Patonga_1min_DPIE" = 'Patonga',          
    "Newcastle_east_1min_PA"  = 'Newcastle East',                
    "PortKembla_BOM_1min_2004" = 'Port Kembla',                
    "PortKembla_1min_BOM" = 'Port Kembla',                   
    "PortKembla_GrainTerminal_1min_PA"  = 'Port Kembla Grain Terminal',      
    "Portland_BOM_1min_2004"  = 'Portland',                
    "Hillarys_BOM_1min_2004" = 'Hillarys',                
    "Hillarys_BOM_1min_2005" = 'Hillarys',                
    "Hillarys_BOM_1min" = 'Hillarys',                     
    "Freemantle_WADoT_5min_2004" = 'Fremantle Boat Harbour',            
    "FreemantleAlternate_WADoT_5min_2004" = 'Fremantle Boat Harbour',   
    "FremantleInnerHarbour_5min" = 'Fremantle Inner Harbour',            
    "ManglesBay_WADoT_5min_2004" = 'Mangles Bay',            
    "BarrackStreet_WADoT_5min_2004" = 'Barrack Street',         
    "MandurahFishermansJetty_WADoT_5min_2004" = "Mandurah Fishermans Jetty",
    "PeelInlet_WADoT_5min_2004" = 'Peel Inlet',             
    "CapeBouvard_WADoT_5min_2004" = 'Cape Bouvard',           
    "Canarvon_WADoT_5min_2004" = 'Carnarvon',              
    "Caddup_WADoT_5min_2004" = 'Caddup',                
    "Harvey_WADoT_5min_2004" = 'Harvey',                
    "Bunbury_WADoT_5min_2004" = 'Bunbury',                
    "BunburyInner_WADoT_5min_2004" = 'Bunbury Inner Harbour',           
    "BusseltonPortGeographe_WADoT_5min_2004" = 'Port Geographe', 
    "Lancelin_WADoT_5min_2004" = 'Lancelin',               
    "JurianBay_WADoT_5min_2004" = 'Jurien Bay',             
    "Geraldton_WADoT_5min_2004" = 'Geraldton',             
    "GeraldtonAlternate_WADoT_5min_2004" = 'Geraldton',     
    "CocosIsland_BOM_1min_2004" = 'Cocos Island',             
    "CNCAR02_01_2004_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2004a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2004b_WADOT_5min" = 'King Bay',           
    "EXEXM02_01_2004_WADOT_5min" = 'Exmouth',            
    "NWONS01_02_2004_WADOT_5min" = 'Onslow Beadon Creek',
    "PWLAM01_03_2004_WADOT_5min" = 'Cape Lambert',            
    "CNCAR02_01_2005_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2005a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2005b_WADOT_5min" = 'King Bay',           
    "EXEXM02_01_2005_WADOT_5min" = 'Exmouth',            
    "NWONS01_02_2005_WADOT_5min" = 'Onslow Beadon Creek',            
    "PWLAM01_03_2005_WADOT_5min" = 'Cape Lambert',            
    "BUBNY01_01_2006_WADOT_5min" = 'Bunbury Inner Harbour',            
    "CNCAR02_01_2006_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2006a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2006b_WADOT_5min" = 'King Bay',           
    "EXEXM02_01_2006_WADOT_5min" = 'Exmouth',            
    "GNGER02_02_2006a_WADOT_5min" = 'Geraldton',           
    "GNGER02_02_2006b_WADOT_5min" = 'Geraldton',           
    "NWONS01_02_2006_WADOT_5min"  = 'Onslow Beadon Creek',           
    "PWLAM01_03_2006_WADOT_5min" = 'Cape Lambert',            
    "BSGEO01_01_2007_WADOT_5min" = 'Port Geographe',            
    "CNCAR02_01_2007_WADOT_5min" = 'Carnarvon',            
    "DPKBY01_03_2007a_WADOT_5min" = 'King Bay',           
    "DPKBY01_03_2007b_WADOT_5min" = 'King Bay',           
    "ESESP03_01_2007_WADOT_5min" = 'Esperance',            
    "EXEXM02_01_2007_WADOT_5min" = 'Exmouth',            
    "GNGER02_02_2007a_WADOT_5min" = 'Geraldton',           
    "GNGER02_02_2007b_WADOT_5min" = 'Geraldton',           
    "NWONS01_02_2007_WADOT_5min" = 'Onslow Beadon Creek',            
    "PWLAM01_03_2007_WADOT_5min" = 'Cape Lambert',
    "NOTHING" = 'Deliberate blank entry'
)

# Get key information in the output df
output_df_provider = unlist(sapply(output_df$sites, match_data_source_to_name, USE.NAMES=FALSE))
output_df_nice_name = rep(NA, nrow(output_df))
for(i in 1:nrow(output_df)) output_df_nice_name[i] = match_site_to_station_name[[output_df$sites[i]]]
combined_df = cbind(output_df, data.frame('site_name'=output_df_nice_name, 'provided_by'=output_df_provider))


# Get the 'nice name' aligned with "events_only"
events_only_nice_name = rep(NA, length(events_only))
for(i in 1:length(events_only)) events_only_nice_name[i] = match_site_to_station_name[[sites_only[i]]]

# New table with a single row per unique site
combined_df_unique_site = combined_df[match(unique(combined_df$site_name), combined_df$site_name),]
# FIX: The previous operation missed the fact that multiple providers give us Fort Dension
k0 = grep("Sydney Fort Denison", combined_df_unique_site$site_name)
combined_df_unique_site$provided_by[k0]='W2018,PANSW'

# Create some codes to describe each event
event_in_order = c('chile1960', 'sumatra2004', 'sumatra2005', 'java2006',
'solomon2007', 'sumatra2007', 'puysegur2009', 'chile2010', 'tohoku2011',
'southamerica2014', 'southamerica2015', 'newhebrides2021', 'kermadec2021',
'sandwich2021')
sites_event_string_data = vector(mode='list', length=nrow(combined_df_unique_site))
for(i in 1:nrow(combined_df_unique_site)){
    sites_event_string_data[[i]] = events_only[which(combined_df_unique_site$site_name[i] == events_only_nice_name)]
    sites_event_string_data[[i]] = sort(unlist(lapply(sites_event_string_data[[i]], function(x) match(x, event_in_order))))
    sites_event_string_data[[i]] = paste0(sites_event_string_data[[i]], collapse=",")
}
sites_event_string_data = unlist(sites_event_string_data)

combined_df_unique_site = cbind(combined_df_unique_site, data.frame(events = sites_event_string_data))
combined_df_unique_site = combined_df_unique_site[order(combined_df_unique_site$site_name),]

# Output in latex format to go in paper
df_for_paper = combined_df_unique_site[c('site_name', 'lon', 'lat', 'events', 'provided_by')]
rownames(df_for_paper) = NULL
library(xtable)
print(xtable(df_for_paper, digits=5))

