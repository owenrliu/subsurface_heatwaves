# Check shear
# Explore velocity shear

library(tidyverse)
library(terra)
library(sf)
library(here)
library(rnaturalearth)

theme_set(theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA)))

# import ocean heat content for the 100-250m depth layer
# one big baby
fls <- list.files(here('data','mom6'),full.names = T) |> str_subset("shear_")
dat <- map(fls,rast) |> rast() |> 
  rotate()

# Make the coastline polygon for measuring distance
coast <- ne_states(country='United States of America',returnclass = 'sf') |> 
  filter(name %in% c('California','Oregon','Washington','Nevada',"Idaho"))
coastcn <- ne_countries(country="Canada",scale='medium',returnclass='sf')
coastmx <- ne_countries(country="Mexico",scale='medium',returnclass='sf')
coast <- st_union(coast,coastcn) |> 
  st_union(coastmx) |> 
  summarise()

# general bounding box
minlat <- 30
maxlat <- 49
minlon <- 230
maxlon <- 244

bb <- st_bbox(c(xmin = minlon, xmax = maxlon, ymin = minlat, ymax = maxlat),
              crs = st_crs(4326))
coast <- st_crop(coast,bb)
ggplot(coast)+geom_sf()
coastv <- vect(coast)

# iso200 <- vect(here('data','200m_isobath_UTM.shp')) |> project(crs(coastv)) |> 
#   as.lines()

# 150km coastline buffer
rt <- dat[[1]]
# rasterize the coastline!
coast_r <- coastv |> 
  # rotate(left=F) |> 
  project(rt) |> 
  # cells on land = 1, else NA
  rasterize(rt,field=1)

# measure the distance of each cell to the continent
dc <- distance(coast_r)/1e3

# Make the raster masks by removing all cells >75 and >150km from coast, respectively
mask150 <- clamp(dc,1,150,values=F) 
plot(mask150)

# mask everything
datm <- mask(dat,mask150)

names(datm) <- paste0("t", seq_len(nlyr(datm)))

bin  <- 0.5                                 # band width in degrees latitude
lat  <- init(datm[[1]], "y")
zone <- floor(lat / bin) * bin + bin / 2
names(zone) <- "lat"

long <- function(r, nm) {
  zonal(r, zone, "mean", na.rm = TRUE) |>
    pivot_longer(-lat, names_to = "layer", values_to = nm)
}

dates <- tibble(layer = names(datm), date = make_date(year=1993,month=1,day=16)+months(0:(nlyr(datm)-1)))

ts <- long(datm, "shear") |>
  left_join(dates, by = "layer") |>
  arrange(lat, date) |> 
  group_by(lat) |> 
  mutate(tidx=row_number()) |> ungroup()

ggplot(ts, aes(tidx, lat, fill = shear)) +
  geom_raster() +
  scale_fill_viridis_c(name = "per s")
