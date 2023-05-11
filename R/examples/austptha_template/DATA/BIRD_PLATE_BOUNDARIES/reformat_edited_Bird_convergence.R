#
# Having edited a shapefile containing Bird's (2003) plate model, we reformat it here.
# This facilitates combining it with the plate motion rates put together by
# Jonathan Griffin for eastern Indonesia and other new source zones, for this project. 
#

library(rptha) # provides readOGR, writeOGR

ptha_conv = readOGR('PTHA_convergence', layer='PTHA_convergence')


local_df = data.frame(
    Long1 = ptha_conv$Long1 + 360*(ptha_conv$Long1 < -40), Lat1 = ptha_conv$Lat1, 
    Long2 = ptha_conv$Long2 + 360*(ptha_conv$Long1 < -40), Lat2 = ptha_conv$Lat2, 
    RL_vel = ptha_conv$RL_vel, Div_vel = ptha_conv$Div_vel, 
    name = rep("", length(ptha_conv)), class = rep("", length(ptha_conv)), 
    Vel_L2R=ptha_conv$Vel_L2R, Azi_Vel=ptha_conv$Azi_Vel,
    collator=rep('Bird2003_subset', length(ptha_conv)))

write.csv(local_df, 'bird_traces_table.csv', row.names=FALSE)


#
# Bird's model does not put convergence directly on Arakan (the northern extension of Sunda Arc),
# but there is evidence that this area is active (Socquet et al., 2006; Cummins, 2007). The global
# strain rate model (if integrated over our source zones in the downdip direction) suggests broadly
# trench perpendicular convergence in this region, although this model allows diffuse deformation over
# a large area and it's possible more of that deformation should be taken up by Arakan (hard to say).
#
arakan_geo = readOGR('Arakan_geometry', layer='Arakan_geometry')
coords = arakan_geo@lines[[1]]@Lines[[1]]@coords

n=50
Longs = approx(coords[,1], n=n+1)$y
Lats  = approx(coords[,2], n=n+1)$y

# Use the Socquet et al (2006) value for Arakan.
library(geosphere)
Azimuth = (bearing( cbind(Longs[1:n], Lats[1:n]), cbind(Longs[2:(n+1)], Lats[2:(n+1)]) ) - 90)%%360
Speed = 23 + Azimuth*0 # mm/year

RL_vel = 0 * Speed
Div_vel = -Speed

new_df = data.frame(
    Long1 = Longs[1:n], Lat1 = Lats[1:n], Long2 = Longs[2:(n+1)], Lat2 = Lats[2:(n+1)],
    RL_vel = RL_vel,
    Div_vel = Div_vel,
    name = rep("Arakan", n),
    class = rep("", n), 
    Vel_L2R = Speed,
    Azi_Vel = Azimuth,
    collator=rep("GD", n))

write.csv(new_df, 'Arakan_traces_table.csv', row.names=FALSE)


#
# Seram-south should have 1mm/year convergence assigned (based on discussion with JG).
# It is currently missing from JG's inputs.
#
seram_south_geo = readOGR('Seramsouth_geometry', layer='Seramsouth_geometry')
coords = seram_south_geo@lines[[1]]@Lines[[1]]@coords

n=50
Longs = approx(coords[,1], n=n+1)$y
Lats  = approx(coords[,2], n=n+1)$y

# Use the Socquet et al (2006) value for Arakan.
library(geosphere)
Azimuth = (bearing( cbind(Longs[1:n], Lats[1:n]), cbind(Longs[2:(n+1)], Lats[2:(n+1)]) ) - 90)%%360
Speed = 1 + Azimuth*0 # mm/year

RL_vel = 0 * Speed
Div_vel = -Speed

new_df = data.frame(
    Long1 = Longs[1:n], Lat1 = Lats[1:n], Long2 = Longs[2:(n+1)], Lat2 = Lats[2:(n+1)],
    RL_vel = RL_vel,
    Div_vel = Div_vel,
    name = rep("seramsouth", n),
    class = rep("", n), 
    Vel_L2R = Speed,
    Azi_Vel = Azimuth,
    collator=rep("GD", n))

write.csv(new_df, 'Seramsouth_traces_table.csv', row.names=FALSE)

