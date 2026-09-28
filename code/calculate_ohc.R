# calculate ocean heat content
library(here)
# library(ncdf4)
library(tidyverse)
library(terra)
library(tictoc)

# url <- "https://compute.earthmover.io/v1/services/dap2/NOAA-PMEL/cefi-nep-hindcast-monthly/main/regrid/main/opendap"
# ds <- nc_open(url)
# 
# # variable names as a table of definitions
# vardefs <- map(names(ds$var),\(x){
#   ncatt_get(ds,x) |> 
#     as_tibble() |> 
#     mutate(varname=x)
# }) |> bind_rows()
# 
# # dimensions as a list of indexed lookup tables
# dims <- imap(ds$dim, \(x,n){
#   tibble(value=c(x$vals),name=n,units=x$units) |> 
#     mutate(index=row_number())
# }) |> set_names(names(ds$dim))
# 
# # times
# timetbl <- dims$time |> 
#   mutate(date=as.Date(value,origin="1993-01-01")) |> 
#   mutate(month=month(date),year=year(date))
# 
# # depths
# dvals <- dims$z_l
# 
# nc_close(ds)

# files with downloaded 3D temperature data
fl <- list.files(here('data','thetao'),full.names = T)
ncdims <- read_rds(here('data','thetao_meta.rds'))
dvals <- ncdims$z_l |> 
  mutate(dbot=value+(lead(value)-value)/2) |> 
  mutate(dbot=coalesce(dbot,6500)) |>
  # layer thickness
  mutate(dz=diff(c(0,dbot)))

# # define a study area bounding box
# latmin <- 32
# latmax <- 48.5
# lonmin <- 230
# lonmax <- 244
# 
# # indices
# lons <- dims$lon |> filter(value >= lonmin, value <= lonmax)
# lats <- dims$lat |> filter(value >= latmin, value <= latmax)
# # subsets of lats/lons
# lon_sub <- pull(lons,value)
# lat_sub <- pull(lats,value)
# # indexes
# lon_idx <-  pull(lons,index)
# lat_idx <- pull(lats,index)
# # counts
# lon_count <- length(lon_sub)
# lat_count=length(lat_sub)
# 
# # a raster extent covering this box, for later plotting
# rext <- ext(min(lon_sub), max(lon_sub), min(lat_sub), max(lat_sub))
# 
# # this helper function makes a raster from a matrix
# # by shaping/dimensioning it correctly
# maker <- function(m){
#   t(m)[rev(seq_len(dim(m)[2])),] |> 
#     rast(
#       extent = rext,
#       crs = "EPSG:4326")
# }

make_ohc <- function(tprof,rho=1026,cp=3990,minz=0,maxz=250){
  depth_range <- c(minz,maxz)
  depth_idx <- which(dvals$value >= depth_range[1] & dvals$value <= depth_range[2])
  dz <- dvals$dz[depth_idx]
  tsub <- tprof[depth_idx]
  # tprof[is.na(tprof)] <- 0
  # integrate across temperature profile
  auc <- sum(tsub * dz, na.rm=FALSE)
  rho * cp * auc / 1e9 # multiply by rho * cp, then convert to gigajoules
}

test <- sds(fl[1])
ohcts <- map(test,\(x) app(x,make_ohc,minz=0,maxz=250))

# Calculate this across all years
# 0-50m
ohcl_z_0_50 <- map(fl,\(x){
  tic(x)
  dat <- sds(x)
  ml <- map(dat,\(d) app(d,make_ohc,minz=0,maxz=50))
  toc()
  rast(ml)
}) |> rast()

# 50-100m
ohcl_z_50_100 <- map(fl,\(x){
  tic(x)
  dat <- sds(x)
  ml <- map(dat,\(d) app(d,make_ohc,minz=50,maxz=100))
  toc()
  rast(ml)
}) |> rast()
# 100-250m
ohcl_z_100_250 <- map(fl,\(x){
  tic(x)
  dat <- sds(x)
  ml <- map(dat,\(d) app(d,make_ohc,minz=100,maxz=250))
  toc()
  rast(ml)
}) |> rast()

# 250-500m
ohcl_z_250_500 <- map(fl,\(x){
  tic(x)
  dat <- sds(x)
  ml <- map(dat,\(d) app(d,make_ohc,minz=250,maxz=500))
  toc()
  rast(ml)
}) |> rast()

# 500-1000m
ohcl_z_500_1000 <- map(fl,\(x){
  tic(x)
  dat <- sds(x)
  ml <- map(dat,\(d) app(d,make_ohc,minz=500,maxz=1000))
  toc()
  rast(ml)
}) |> rast()


# save with correct times/dates
timetbl <- ncdims$time |> 
  mutate(date=as.Date(value,origin="1993-01-01")) |> 
  mutate(month=month(date),year=year(date))
time(ohcl_z_0_50) <- time(ohcl_z_50_100) <- time(ohcl_z_100_250) <- time(ohcl_z_250_500) <- time(ohcl_z_500_1000) <- timetbl$date

writeRaster(ohcl_z_0_50,here('output','ocean_heat_content_0_50_ts.tif'),overwrite=T)
writeRaster(ohcl_z_50_100,here('output','ocean_heat_content_50_100_ts.tif'),overwrite=T)
writeRaster(ohcl_z_100_250,here('output','ocean_heat_content_100_250_ts.tif'),overwrite=T)
# writeRaster(ohcl_z_0_250,here('output','ocean_heat_content_0_250_ts.tif'),overwrite=T)
writeRaster(ohcl_z_250_500,here('output','ocean_heat_content_250_500_ts.tif'),overwrite=T)
writeRaster(ohcl_z_500_1000,here('output','ocean_heat_content_500_1000_ts.tif'),overwrite=T)
