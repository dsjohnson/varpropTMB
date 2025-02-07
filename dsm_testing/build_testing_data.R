
setwd("~/research/projects/r_packages/varpropTMB/dsm_testing")

library(dsm)
library(ggplot2)

# plotting options
gg.opts <- theme(panel.grid.major=element_blank(),
                 panel.grid.minor=element_blank(),
                 panel.background=element_blank())
load("mexdolphins-extra.rda")
data(mexdolphins)
library(sf)
library(plyr)

# tell R that the survey.area object is currently in lat/long
sp::proj4string(survey.area) <- sp::CRS("+proj=longlat +datum=WGS84")

predsf <- st_as_sf(pred.polys)

area.sf <- st_as_sf(survey.area)
st_crs(area.sf) <- "WGS84"
area.sf.proj <- st_transform(area.sf, crs = st_crs(predsf))

# Convert preddata to a spatial object
preddata_sf <- st_as_sf(preddata, coords=c("x", "y"))
st_crs(preddata_sf) <- st_crs(area.sf.proj)
# Perform the intersection
preddata_sf <- st_intersection(preddata_sf, area.sf.proj)
coords_preddata <- data.frame(st_coordinates(preddata_sf))

preddata_sf$x <- coords_preddata$X
preddata_sf$y <- coords_preddata$Y 
# proj 4 string
# using http://spatialreference.org/ref/esri/north-america-lambert-conformal-conic/
lcc_proj4 <- sp::CRS("+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")

# project using LCC
survey.area <- sp::spTransform(survey.area, CRSobj=lcc_proj4)

# simplify the object
survey.area <- data.frame(survey.area@polygons[[1]]@Polygons[[1]]@coords)
names(survey.area) <- c("x", "y")
segdata_sf <- st_as_sf(segdata, coords = c("x","y"))
st_crs(segdata_sf) <- st_crs(area.sf.proj)
# study area outline and segment centres
plot_segments <- ggplot() +
  geom_sf(data = area.sf.proj, fill="lightblue", color = "blue", linewidth=.1) +
  geom_sf(data=segdata_sf, fill=NA, color="black", linewidth=.3) +
  labs(title="1996 SE Fisheries Science Center Gulf of Mexico cruise",
       subtitle = "Points are segment centres") +
  scale_fill_viridis_c(option = "magma", guide = "none")
plot_segments

predsf <- st_as_sf(pred.polys)
# plot as projected
plot(st_geometry(predsf), axes=TRUE)

prediction_grid <- st_make_grid(area.sf.proj, cellsize = c(9000,9000))
prediction_grid_sf <- st_sf(geometry = prediction_grid)
cropped_grid <- st_join(prediction_grid_sf, preddata_sf, join = st_nearest_feature)
cropped_grid <- st_intersection(cropped_grid, area.sf.proj)

depth <- ggplot() +
  geom_sf(data=cropped_grid, aes(fill=depth), color=NA) +
  labs(title = "Spotted dolphins, Gulf of Mexico",
       subtitle = "Depth in meters, size of detected dolphin groups") +
  xlab("Longitude") + ylab("Latitude") +
  geom_point(aes(x, y, size=size), data=distdata, colour="red",alpha=I(0.5)) +
  scale_fill_viridis_c(option = "viridis", direction = 1)
depth

library(Distance)
detfc.hr.null <- ds(distdata, max(distdata$distance), key="hr", adjustment=NULL)

detfc.hr.beau<-ds(distdata, max(distdata$distance), formula=~as.factor(beaufort),
                  key="hr", adjustment=NULL)



dsm.xy <- dsm(count~s(x,y), detfc.hr.null, segdata, obsdata, method="REML")

summary(dsm.xy)

vis.gam(dsm.xy, plot.type="contour", view=c("x","y"), asp=1, type="response", contour.col="black", n.grid=100)

dsm.xy.depth <- dsm(count~s(x,y,k=10) + s(depth,k=20), detfc.hr.null, segdata, obsdata, method="REML")
summary(dsm.xy.depth)

dsm.est.xy <- dsm(abundance.est~s(x,y), detfc.hr.beau, segdata, obsdata, method="REML")

dsm.xy.tweedie <- dsm(count~s(x,y), detfc.hr.null, segdata, obsdata, family=tw(), method="REML")
summary(dsm.xy.tweedie)

dsm.xy.pred <- predict(dsm.xy, preddata, preddata$area)

prediction_grid <- st_make_grid(area.sf.proj, cellsize = c(9000,9000))
prediction_grid_sf <- st_sf(geometry = prediction_grid)

preddata_sf$Prediction_xy <- dsm.xy.pred

cropped_grid <- st_join(prediction_grid_sf, preddata_sf, join = st_nearest_feature)
cropped_grid <- st_intersection(cropped_grid, area.sf.proj)

dsm.xy.depth.pred <- predict(dsm.xy.depth, preddata, preddata$area)
preddata_sf$Prediction_xy_depth <- dsm.xy.depth.pred
cropped_grid <- st_join(prediction_grid_sf, preddata_sf, join = st_nearest_feature)
cropped_grid <- st_intersection(cropped_grid, area.sf.proj)


save(distdata, segdata, obsdata, preddata, segdata_sf, area.sf, preddata_sf, detfc.hr.null, dsm.xy.depth, file="dsm_test_data.RData")


