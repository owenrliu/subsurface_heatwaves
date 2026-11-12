# Explore OHC gradient

library(tidyverse)
library(terra)
library(sf)
library(here)
library(rnaturalearth)

theme_set(theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA)))

# import ocean heat content for the 100-250m depth layer
dat <- rast(here('output','ocean_heat_content_100_250_ts.tif')) |> rotate()

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

iso200 <- vect(here('data','200m_isobath_UTM.shp')) |> project(crs(coastv)) |> 
  as.lines()

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

# Sobel kernels (weights for the focal calc, so result is a per-cell derivative)
kx <- matrix(c(-1, 0, 1,
               -2, 0, 2,
               -1, 0, 1), 3, byrow = TRUE) / 8
ky <- t(kx) * -1                  # positive = increasing northward

# dTdx_cell <- focal(dat, kx, fun = "sum", na.policy = "omit", fillvalue = NA)
# dTdy_cell <- focal(dat, ky, fun = "sum", na.policy = "omit", fillvalue = NA)

dTdx_cell <- focal(datm, kx, fun = "sum", fillvalue = NA)
dTdy_cell <- focal(datm, ky, fun = "sum",  fillvalue = NA)


# Cell size in km, varying with latitude
lat   <- init(datm, "y")
dx_km <- 111.32 * cos(lat * pi / 180) * res(datm)[1]
dy_km <- 110.57 * res(datm)[2]

gx <- dTdx_cell / dx_km           # OHC/km, eastward
gy <- dTdy_cell / dy_km           # OHC/km, northward

grad_mag <- sqrt(gx^2 + gy^2)     # gradient magnitude, OHC/km
grad_dir <- atan2(gy, gx) * 180 / pi   # direction of gradient, degrees from east

plot(grad_mag[[1]], main = "OHC gradient magnitude (GJ/km)")
plot(grad_dir[[1]], main = "OHC gradient direction (GJ/km)")

## mapping clearly
thr <- global(grad_mag, quantile, probs = 0.5, na.rm = TRUE)[1, 1]
# only plot strong magnitudes
dir_strong <- ifel(grad_mag >= thr, grad_dir %% 360, NA)

plot(dir_strong[[254]], range = c(0, 360),
     col = hcl(h = seq(0, 360, length.out = 360), c = 70, l = 65),
     main = "Direction of OHC gradient, strongest 50% of gradients")

plot_grad_arrows <- function(magr,agg=4,tidx){
    
    # Fraction of fine cells in each block that are valid
    valid  <- !is.na(magr)
    frac   <- terra::aggregate(valid, fact=agg, "mean")
    
    gxa <- terra::aggregate(gx, agg, "mean", na.rm = TRUE)
    gya <- terra::aggregate(gy, agg, "mean", na.rm = TRUE)
    
    # Drop blocks with less than, say, 70% valid cells
    # gxa <- ifel(frac >= 0.6, gxa, NA)
    # gya <- ifel(frac >= 0.6, gya, NA)
    
    v <- as.data.frame(c(gxa[[tidx]], gya[[tidx]]), xy = TRUE, na.rm = TRUE)
    names(v) <- c("x", "y", "gx", "gy")
    v$mag <- sqrt(v$gx^2 + v$gy^2)
    v <- v[v$mag >= thr, ]
    
    L    <- agg * res(datm)[2] * 0.8
    v$x2 <- v$x + (v$gx / v$mag) * L / cos(v$y * pi / 180)
    v$y2 <- v$y + (v$gy / v$mag) * L
    
    # Keep only arrows whose tail AND head fall on valid water
    ok_tail <- !is.na(extract(magr[[tidx]], v[, c("x", "y")])[, 2])
    ok_head <- !is.na(extract(magr[[tidx]], v[, c("x2", "y2")])[, 2])
    v <- v[ok_tail & ok_head, ]
    
    timelab <- time(dat)[tidx] |> as.character()
    
    plot(magr[[tidx]],col=map.pal("water",100),colNA="gray90",main=timelab)
    arrows(v$x, v$y, v$x2, v$y2, length = 0.05, lwd = 1.5,col="black")
    plot(coast,add=T,border='black',lwd=1)
    # plot(iso200,add=T,col='black')
}

plot_grad_arrows(magr=grad_mag,agg=3,tidx=240)
plot_grad_arrows(magr=grad_mag,agg=3,tidx=250)
plot_grad_arrows(magr=grad_mag,agg=3,tidx=260)

# rose diagram

d <- na.omit(data.frame(dir = as.vector(values(grad_dir %% 360)),
                        mag = as.vector(values(grad_mag))))
d <- d[d$mag >= thr, ]

ggplot(d, aes(dir, weight = mag)) +
  geom_histogram(binwidth = 15, boundary = 0, fill = "grey40") +
  coord_polar(start = pi / 2, direction = -1) +   # 0° at east, counterclockwise
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 315, 45))

# Time series

names(grad_mag) <- names(grad_dir) <- paste0("t", seq_len(nlyr(grad_mag)))

gx <- grad_mag * cos(grad_dir * pi / 180)   # grad_dir: degrees from east, as before
gy <- grad_mag * sin(grad_dir * pi / 180)

bin  <- 0.5                                 # band width in degrees latitude
lat  <- init(grad_mag[[1]], "y")
zone <- floor(lat / bin) * bin + bin / 2
names(zone) <- "lat"

long <- function(r, nm) {
  zonal(r, zone, "mean", na.rm = TRUE) |>
    pivot_longer(-lat, names_to = "layer", values_to = nm)
}

dates <- tibble(layer = names(grad_mag), date = time(grad_mag))

ts <- long(grad_mag, "mag") |>
  left_join(long(gx, "gx"), by = c("lat", "layer")) |>
  left_join(long(gy, "gy"), by = c("lat", "layer")) |>
  left_join(dates, by = "layer") |>
  mutate(
    dir_math    = (atan2(gy, gx) * 180 / pi) %% 360,   # 0 = east, counterclockwise
    dir_bearing = (90 - dir_math) %% 360,              # compass bearing of warming direction
    consistency = sqrt(gx^2 + gy^2) / mag              # 0 = directions cancel, 1 = all aligned
  ) |>
  arrange(lat, date) |> 
  group_by(lat) |> 
  mutate(tidx=row_number()) |> ungroup()

ggplot(ts, aes(tidx, lat, fill = mag)) +
  geom_raster() +
  scale_fill_viridis_c(name = "GJ/km")

ggplot(filter(ts, consistency >= 0.3), aes(tidx, lat, fill = dir_bearing)) +
  geom_raster() +
  scale_fill_gradientn(colours = hcl(seq(0, 360, length.out = 360), 70, 65),
                       limits = c(0, 360), name = "Bearing")
# climatology
clim <- ts |>
  mutate(month = factor(as.integer(format(date, "%m")), 1:12, month.abb)) |>
  group_by(lat, month) |>
  summarise(mag = mean(mag, na.rm = TRUE),
            gx  = mean(gx,  na.rm = TRUE),
            gy  = mean(gy,  na.rm = TRUE),
            n   = sum(!is.na(mag)), .groups = "drop") |>
  mutate(dir_math    = (atan2(gy, gx) * 180 / pi) %% 360,
         dir_bearing = (90 - dir_math) %% 360,
         consistency = sqrt(gx^2 + gy^2) / mag) |>
  group_by(lat) |>
  mutate(mag_anom = mag - mean(mag, na.rm = TRUE)) |>   # seasonal departure within each band
  ungroup()
# Magnitude climatology
ggplot(clim, aes(month, lat, fill = mag)) +
  geom_tile() + scale_fill_viridis_c(name = "GJ/km")

# Seasonal cycle with the latitudinal mean removed (shows when each band peaks)
ggplot(clim, aes(month, lat, fill = mag_anom)) +
  geom_tile() +
  scale_fill_gradient2(name = "Anomaly\n(GJ/km)", midpoint = 0)

# Direction, masked where the monthly vectors cancel out
ggplot(filter(clim, consistency >= 0.3), aes(month, lat, fill = dir_bearing)) +
  geom_tile() +
  scale_fill_gradientn(colours = hcl(seq(0, 360, length.out = 360), 70, 65),
                       limits = c(0, 360), name = "Bearing")

# time series around Cape Blanco
ts |> 
  filter(lat>41.5,lat<43.5) |> 
  ggplot(aes(date,dir_bearing,color=factor(lat)))+
  facet_wrap(~factor(lat))+
  geom_line()

clim |> 
  filter(lat>41.5,lat<43.5) |> 
  mutate(lat=as.factor(lat)) |> 
  ggplot(aes(month,dir_bearing,color=lat,group=lat))+
  # facet_wrap(~factor(lat))+
  geom_line()+
  labs(x="Month",y="Gradient Direction")

clim |> 
  filter(lat>41.5,lat<43.5) |> 
  mutate(lat=as.factor(lat)) |> 
  ggplot(aes(month,consistency,color=lat,group=lat))+
  # facet_wrap(~factor(lat))+
  geom_line()+
  labs(x="Month",y="Gradient Consistency")

clim |> 
  filter(lat>41.5,lat<43.5) |> 
  mutate(lat=as.factor(lat)) |> 
  ggplot(aes(month,mag,color=lat,group=lat))+
  # facet_wrap(~factor(lat))+
  geom_line()+
  labs(x="Month",y="Gradient Magnitude")

# remove clim to calc anoms
climj <- clim |> 
  select(lat,month,mag,dir_bearing) |> 
  rename(magclim=mag,dirclim=dir_bearing)
tsanom <- ts |> 
  mutate(month = factor(as.integer(format(date, "%m")), 1:12, month.abb)) |>
  left_join(climj,by = join_by(lat, month)) |>
  mutate(magdiff=(mag-magclim),bearingdiff=dir_bearing-dirclim)

tsanom |> 
  filter(lat>41.5,lat<43.5) |> 
  ggplot(aes(date,magdiff,color=factor(lat),alpha=consistency))+
  facet_wrap(~factor(lat))+
  geom_hline(yintercept=0,linetype=2)+
  geom_line()

tsanom |> 
  filter(lat>41.5,lat<43.5) |> 
  ggplot(aes(date,bearingdiff,color=factor(lat),alpha=consistency))+
  facet_wrap(~factor(lat))+
  geom_point()

tsanom |> 
  filter(consistency>0.3) |> 
  mutate(year=year(date)) |> 
  filter(lat>41.5,lat<43.5) |> 
  mutate(lat=as.factor(lat)) |> 
  ggplot(aes(month,dir_bearing,group=year))+
  facet_wrap(~factor(lat))+
  geom_point()+
  labs(x="Month",y="Gradient Direction")+
  theme(axis.text.x=element_text(angle=45,hjust=1,vjust=1))

tsanom |> 
  ggplot(aes(bearingdiff,magdiff))+geom_point()

tsanom |> 
  filter(lat>41.5,lat<43.5) |>
  ggplot(aes(bearingdiff))+geom_density()


# calculate dominant direction per band
library(circular)

circ_summary <- function(dir_deg, w, kappa = 8) {
  ok <- is.finite(dir_deg) & is.finite(w) & w > 0
  th <- dir_deg[ok] * pi / 180
  w  <- w[ok]
  if (length(th) < 5) return(tibble(n = length(th)))
  
  res <- function(a) c(C = sum(w * cos(a)), S = sum(w * sin(a))) / sum(w)
  v1 <- res(th); v2 <- res(2 * th)
  
  grid <- seq(0, 2 * pi, length.out = 361)[-361]
  dens <- vapply(grid, \(g) sum(w * exp(kappa * cos(g - th))), numeric(1))
  
  tibble(
    n           = length(th),
    mean_dir    = (atan2(v1["S"], v1["C"]) * 180 / pi) %% 360,
    R           = sqrt(sum(v1^2)),
    mode_dir    = grid[which.max(dens)] * 180 / pi,
    R_axial     = sqrt(sum(v2^2)),
    axis_dir    = ((atan2(v2["S"], v2["C"]) / 2) * 180 / pi) %% 180,   # orientation, 0-180
    rayleigh_p  = rayleigh.test(circular(th, units = "radians"))$p.value
  )
}

dom <- ts |>
  mutate(month = factor(as.integer(format(date, "%m")), 1:12, month.abb)) |>
  group_by(lat,month) |>
  reframe(circ_summary(dir_math, mag)) |>
  mutate(across(c(mean_dir, mode_dir), \(x) (90 - x) %% 360, .names = "{.col}_brg"))

dom <-dom |> 
  left_join(select(clim,lat,month,mag))

ggplot(dom, aes(lat, mode_dir_brg, alpha = R)) +
  facet_wrap(~month)+
  geom_point() +
  scale_y_continuous(limits = c(0, 360), breaks = seq(0, 360, 90),
                     labels=c("0(N)","90(E)","180(S)","270(W)","360(N)")) +
  # scale_colour_viridis_c() +
  labs(y = "Dominant bearing of OHC (°)", x = "Latitude")+
  theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA))

ggplot(dom, aes(month, mode_dir_brg,group=lat)) +
  facet_wrap(~lat)+
  geom_point(aes(alpha=R)) +
  geom_line(linetype=2)+
  scale_y_continuous(limits = c(0, 360), breaks = seq(0, 360, 90),
                     labels=c("0(N)","90(E)","180(S)","270(W)","360(N)")) +
  # scale_colour_viridis_c() +
  labs(y = "Dominant bearing of OHC (°)", x = "Latitude")+
  theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA))

ggplot(dom, aes(month, mag,group=lat)) +
  facet_wrap(~lat)+
  geom_point(aes(alpha=R)) +
  geom_line(linetype=2)+
  # scale_y_continuous(limits = c(0, 360), breaks = seq(0, 360, 90),
  #                    labels=c("0(N)","90(E)","180(S)","270(W)","360(N)")) +
  # scale_colour_viridis_c() +
  labs(y = "Gradient Magnitude", x = "Latitude")+
  theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA),
                        axis.text.x=element_text(angle=45,hjust=1,vjust=1))

ggplot(dom, aes(x = mode_dir_brg, y = lat, colour = R, size = R)) +
  geom_point(alpha = 0.8) +
  coord_polar(start = 0, direction = 1) +           # 0° at top, clockwise
  scale_x_continuous(limits = c(0, 360),
                     breaks = seq(0, 315, 45),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) +
  scale_y_continuous(limits = c(min(dom$lat) - 2, NA)) +   # leaves a hole at the center
  scale_colour_viridis_c(limits = c(0, 1)) +
  scale_size_continuous(range = c(1, 5), limits = c(0, 1)) +
  guides(colour = guide_legend(), size = guide_legend()) +
  labs(x = NULL, y = "Latitude (°N)",
       colour = "R", size = "R",
       title = "Dominant bearing of gradient by latitude") +
  theme_minimal()

ggplot(dom, aes(x = mode_dir_brg, y = R, colour = lat)) +
  geom_point(size = 2.5) +
  coord_polar(start = 0, direction = 1) +
  scale_x_continuous(limits = c(0, 360), breaks = seq(0, 315, 45),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_colour_viridis_c(name = "Latitude (°N)") +
  labs(x = NULL, y = "R") +
  theme_minimal()

# looking at how consistent vs. multimodal the bearings are
dom |> ggplot(aes(R,R_axial))+geom_point()+geom_abline(slope=1,intercept=0)

dom |>   
  ggplot(aes(x = mode_dir_brg,y=mag, colour = lat)) +
  geom_segment(aes(xend=mode_dir_brg,yend=mag,y=0),arrow = arrow(length=unit(0.1,"cm")))+
  # geom_point(alpha = 0.8) +
  # geom_col()+
  coord_polar(start = 0, direction = 1) +           # 0° at top, clockwise
  scale_x_continuous(limits = c(0, 360),
                     breaks = seq(0, 315, 45),
                     labels = c("N", "NE", "E", "SE", "S", "SW", "W", "NW")) +
  scale_y_continuous(limits = c(min(dom$mag) - 0.002, NA)) +   # leaves a hole at the center
  scale_colour_viridis_c(name = "Latitude (°N)") +
  # scale_size_continuous(range = c(1, 5), limits = c(0, 1)) +
  guides(colour = guide_legend(), size = guide_legend()) +
  labs(x = NULL, y = "Magnitude",
       colour = "Latitude (°N)", size = "Latitude (°N)",
       title = "Gradient magnitude by latitude") +
  theme_minimal()
