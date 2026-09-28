# download daily earthmover MOM6 temperature data
library(ncdf4)
# library(tidync)
library(tidyverse)
library(terra)
library(tictoc)
library(here)

url <- "https://compute.earthmover.io/v1/services/dap2/NOAA-PMEL/cefi-nep-hindcast-monthly/main/regrid/main/opendap"
ds <- nc_open(url)

# variable names as a table of definitions
vardefs <- map(names(ds$var),\(x){
  ncatt_get(ds,x) |> 
    as_tibble() |> 
    mutate(varname=x)
}) |> bind_rows()

# dimensions as a list of indexed lookup tables
dims <- imap(ds$dim, \(x,n){
  tibble(value=c(x$vals),name=n,units=x$units) |> 
    mutate(index=row_number())
}) |> set_names(names(ds$dim))
write_rds(dims,here('data','thetao_meta.rds'))

timetbl <- dims$time |> 
  mutate(date=as.Date(value,origin="1993-01-01")) |> 
  mutate(month=month(date),year=year(date))

# define a study area bounding box
latmin <- 32
latmax <- 48.5
lonmin <- 230
lonmax <- 244

# indices
lons <- dims$lon |> filter(value >= lonmin, value <= lonmax)
lats <- dims$lat |> filter(value >= latmin, value <= latmax)
# subsets of lats/lons
lon_sub <- pull(lons,value)
lat_sub <- pull(lats,value)
# indexes
lon_idx <-  pull(lons,index)
lat_idx <- pull(lats,index)
# counts
lon_count <- length(lon_sub)
lat_count=length(lat_sub)

# a raster extent covering this box, for later plotting
rext <- ext(min(lon_sub), max(lon_sub), min(lat_sub), max(lat_sub))

# convenience function to reshape a time-slice and make it into a raster
maker <- function(m){
  t(m)[rev(seq_len(dim(m)[2])),] |> 
    rast(
      extent = rext,
      crs = "EPSG:4326")
}

# download function
download_nc_yr <- function(ncurl=url,yr,min_zind=1,max_zind=36){ # 0-1000m default
  # find time indices
  ti <- timetbl |> filter(year==yr) |> pull(index)
  dt <- timetbl$date[ti]
  tstart <- min(ti)
  tcount <- length(ti)
  
  # get temperature
  # tic(paste("Pulling",length(min_zind:max_zind),"temperature slices"))
  nc <- nc_open(ncurl)
  thetao <- ncvar_get(
    nc,
    varid="thetao",
    # xyzt (time is last dimension)
    start=c(min(lon_idx),min(lat_idx),min_zind,tstart),
    count=c(lon_count,lat_count,length(min_zind:max_zind),tcount)
  )
  # toc()
  
  # make raster
  nt <- dim(thetao)[4] # times (should be 12 months in a year)
  nz <- dim(thetao)[3]
  
  ml <- list()
  for(i in 1:nt){
    m <- list()
    for(j in 1:nz){
      m[[j]] <- maker(thetao[,,j,i])
      time(m[[j]]) <- dt[i]
    }
    mr <- rast(m)
    ml[[i]] <- mr
  }

  rthetao <- sds(ml)
  names(rthetao) <- month.name[month(dt)]
  
  fn <- paste0("thetao_",yr,".nc")
  
  writeCDF(rthetao,here('data','thetao',fn),varname="thetao",
           longname="Sea Water Potential Temperature",
           overwrite=T)
  
  nc_close(nc)
}

# download_nc_yr(yr=1997)
# test <- sds(here('data','thetao','thetao_2025.nc'))

for(i in unique(timetbl$year)){
  tic(paste("Downloading",i,"done"))
  download_nc_yr(yr=i)
  toc()
}
