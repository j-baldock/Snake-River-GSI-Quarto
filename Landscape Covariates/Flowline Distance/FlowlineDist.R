
#----------------------------------------------------------------------------------------------------#
# Calculate distances between reporting unit locations and mainstem Snake River section midpoints
#----------------------------------------------------------------------------------------------------#

library(tidyverse)
library(riverdist)
library(mapview)
library(sf)



# load flowline as riverdist network
flowline <- line2network(path = "/Users/jeffbaldock/Library/CloudStorage/GoogleDrive-jbaldock@uwyo.edu/Shared drives/wyo-coop-baldock/UWyoming/Snake River Cutthroat/Data/Spatial/Flowline/FlowlineEditing", layer = "SnakeHeadwaters_flowline_springsclean")
plot(flowline)
zoomtoseg(seg = c(67,3423), rivers = flowline)

# fix topological errors (from raw flowline shp, mouth at segment 65 and vertex 94)
flowline_fixed <- cleanup(flowline)
save(flowline_fixed, file = "Landscape Covariates/Flowline Distance/SnakeHeadwaters_RiverDist_Cleaned.Rdata")
load(file = "Landscape Covariates/Flowline Distance/SnakeHeadwaters_RiverDist_Cleaned.Rdata")
# cleanup_verts(flowline_fixed)
# topologydots(flowline_fixed)

# load sites and convert to shapefile
sites_repu <- sf::read_sf(dsn = "Landscape Covariates", layer = "RepUnit_SpatialLocations")
sites_mixt <- read_csv("Landscape Covariates/Snake_section_midpoints.csv")
sites_mixt1 <- sf::st_as_sf(sites_mixt, coords = c("lon", "lat"), crs = 4326) %>% sf::st_transform(crs = sf::st_crs(sites_repu))
names(sites_mixt1) <- c("repunit", "geometry")

# convert to river points and snap to network
sites_repu2 <- xy2segvert(x = sf::st_coordinates(sites_repu)[,1], y = sf::st_coordinates(sites_repu)[,2], rivers = flowline_fixed)
sites_mixt2 <- xy2segvert(x = sf::st_coordinates(sites_mixt1)[,1], y = sf::st_coordinates(sites_mixt1)[,2], rivers = flowline_fixed)

# inspect
zoomtoseg(seg = c(2077), rivers = flowline_fixed)
points(sf::st_coordinates(sites_repu)[,1], sf::st_coordinates(sites_repu)[,2], pch = 16, col = "red")
riverpoints(seg = sites_repu2$seg, vert = sites_repu2$vert, rivers = flowline_fixed, pch = 16, col = "blue")

# compute distances between reporting units and Snake River section mid-points 
distlist <- list()
for (i in 1:dim(sites_mixt1)[1]) {
  sites <- rbind(sites_repu, sites_mixt1[i,])
  sites2 <- xy2segvert(x = st_coordinates(sites)[,1], y = st_coordinates(sites)[,2], rivers = flowline_fixed)
  disttib <- tibble(section = rep(sites_mixt1$repunit[i], times = 57), repunit = NA, distkm = NA)
  for (j in 1:52) {
    disttib$repunit[j] <- sites$repunit[j]
    distm <- riverdistance(startseg = sites2$seg[53], startvert = sites2$vert[53],
                           endseg = sites2$seg[j], endvert = sites2$vert[j], rivers = flowline_fixed, map = F)
    disttib$distkm[j] <- distm/1000
    print(j)
  }
  distlist[[i]] <- disttib
}
disttib <- do.call(rbind, distlist)
write_csv(disttib, "Landscape Covariates/Flowline Distance/SnakeRiverSections_RepUnits_FlowlineDistance.csv")
disttib <- read_csv("Landscape Covariates/Flowline Distance/SnakeRiverSections_RepUnits_FlowlineDistance.csv")

# plot histogram
jpeg("Landscape Covariates/Flowline Distance/FlowlineDistance_hist.jpg", units = "in", res = 1000, width = 5, height = 5)
hist(disttib$distkm, xlab = "Flowline distance (km)", main = "Distance between reporting units and \nmainstem Snake River section midpoints")
dev.off()
