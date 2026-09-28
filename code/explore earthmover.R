# Lazily load the entire dataset
library(ncdf4)
# library(tidync)
library(tidyverse)
library(terra)
library(tictoc)
url <- "https://compute.earthmover.io/v1/services/dap2/NOAA-PMEL/cefi-nep-hindcast-monthly/main/regrid/main/opendap"
ds <- nc_open(url)
theme_set(theme_minimal()+theme(panel.border = element_rect(color='black',fill=NA)))

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
timetbl <- dims$time |> 
  mutate(date=as.Date(value,origin="1993-01-01")) |> 
  mutate(month=month(date),year=year(date))
# e.g. for all the latitude values, table dims$lat

# define a study area bounding box
latmin <- 25
latmax <- 60
lonmin <- 200
lonmax <- max(dims$lon$value)

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

# grab a slice as a function of a time and depth index
get_slice <- function(ncin,varname,tidx,didx){
  out <- ncvar_get(
    nc=ncin,
    varid=varname,
    # xyzt (time is last dimension)
    start=c(min(lon_idx),min(lat_idx),didx,tidx),
    count=c(lon_count,lat_count,1,1)
  )
  out
}
test <- get_slice(ds,"thetao",1,1)

# raster extent
rext <- ext(min(lon_sub), max(lon_sub), min(lat_sub), max(lat_sub))
# make sure we're filling the raster correctly (lats change fastest, then lons)
# and north to south
testf <- t(test)[rev(seq_len(ncol(test))),]

r <- rast(
  testf,
  extent = rext,
  crs = "EPSG:4326"
)

plot(r)

r_slice <- function(ncin,varname,tidx,didx){
  d <- get_slice(ncin,varname,tidx,didx)
  r <- t(d)[rev(seq_len(ncol(d))),] |> 
    rast(
      extent = rext,
      crs = "EPSG:4326")
  r
}

r_slice(ncin=ds,varname="thetao",tidx=1,didx=21) |> plot()

# try finding isopycnals
# find a density variable
vardefs[grepl("density",vardefs$long_name),]
# rhopot0
rhop <- r_slice(ncin=ds,varname="rhopot0",tidx=1,didx=36)
plot(rhop)

# try to pull all depth slices of a time point
tidx=1;min_zind=1;max_zind=36;ncin <- ds;target_rho=1022
rm(tidx,min_zind,max_zind,ncin)

get_3Dslice <- function(ncurl,varname,tidx,min_zind,max_zind){
  tic("opening url")
  nc <- nc_open(ncurl)
  toc()
  tic("pulling variable")
  out <- ncvar_get(
    nc,
    varid=varname,
    # xyzt (time is last dimension)
    start=c(min(lon_idx),min(lat_idx),min_zind,tidx),
    count=c(lon_count,lat_count,length(min_zind:max_zind),1)
  )
  toc()
  out
}

test <- get_3Dslice(url,varname="thetao",tidx=10,min_zind=1,max_zind=36)

calc_isopyc <- function(ncurl,tidx,target_rho,min_zind=1,max_zind=52,rast=F){
  nc <- nc_open(ncurl)
  tic(paste("Pulling",length(min_zind:max_zind),"depth slices"))
  rhopot <- ncvar_get(
    nc,
    varid="rhopot0",
    # xyzt (time is last dimension)
    start=c(min(lon_idx),min(lat_idx),min_zind,tidx),
    count=c(lon_count,lat_count,length(min_zind:max_zind),1)
  )
  toc()
  
  tic(paste("Pulling",length(min_zind:max_zind),"temperature slices"))
  thetao <- ncvar_get(
    nc,
    varid="thetao",
    # xyzt (time is last dimension)
    start=c(min(lon_idx),min(lat_idx),min_zind,tidx),
    count=c(lon_count,lat_count,length(min_zind:max_zind),1)
  )
  toc()
  
  # find which depths are closest to the target density
  tic("Calculating isopycnal")
  which_depths <- apply(rhopot, c(1, 2), function(profile) {
    vals <- abs(profile - target_rho)
    minval <- min(vals)
    ind <- which.min(vals)
    if(length(ind)==0) ind <- NA
    if(is.na(minval)|minval>0.2) ind <- NA
    ind
  })
  which_depths
  toc()
  
  # extract temperature at those depths
  # Extract temperature (need to convert indices back to (i, j, k) form)
  tic("Extracting temperature along isopycnal")
  d <- matrix(NA, nrow=nrow(rhopot), ncol=ncol(rhopot))
  for (i in 1:nrow(rhopot)) {
    for (j in 1:ncol(rhopot)) {
      k <- unlist(which_depths[i, j])
      d[i, j] <- thetao[i, j, k]
    }
  }
  toc()
  if(rast){
    r <- t(d)[rev(seq_len(ncol(d))),] |> 
      rast(
        extent = rext,
        crs = "EPSG:4326")
    return(r)
  } else return(d)
}

test <- calc_isopyc(ncurl=url,min_zind=1,max_zind=36,tidx = 1,target_rho = 1027)
plot(test)

## Strategy- use the data to try to identify depths (regionally) with the greatest temperature variance
## Then find isopycnals associated with those depths
## Then track heatwaves along those isopycnals

# pick a few random timepoints and grab 3D temperature
set.seed(1234)
tsamp <- sample(dims$time$index,size = 10,replace=F)
# pick 2014?
tsamp <- timetbl |> filter(year==2014) |> pull(index)
t3d <- map(tsamp,get_3Dslice,ncurl=url,varname='thetao',min_zind=1,max_zind=36)
t4d <- simplify2array(t3d)
# for each sample timepoint, calculate mean temperature profile

profile_t <- function(tmat,depths){
  temp_sd <- apply(tmat,3,function(x) sd(as.vector(x),na.rm=T))
  mean_temp <- apply(tmat,3,mean,na.rm=T)
  tibble(depth=depths,mean_temp=mean_temp,sd_temp=temp_sd)
}

depths <- dims$z_l |> arrange(index) |> slice(1:36) |> pull(value)
t3df <- imap(t3d,\(x,y) profile_t(x,depths=depths) |> mutate(samp=y)) |> bind_rows()

p_sd_temp <- t3df |> ggplot(aes(depth,sd_temp,color=factor(samp)))+geom_line()+
  labs(x="Depth",y="SD Across Spatial Domain",color="Month")
p_mean_temp <- t3df |> ggplot(aes(depth,mean_temp,color=factor(samp)))+geom_line()+
  labs(x="Depth",y="Mean Across Spatial Domain",color="Month")

p_mean_temp
p_sd_temp

# Find depth(s) with highest temperature variance
depth_max_var <- t3df |> group_by(samp) |> 
  slice_max(sd_temp) |> 
  ungroup()
depth_max_var |> ggplot(aes(samp,depth))+geom_line()+
  scale_x_continuous(breaks=seq_len(12))+
  labs(y="Depth of Max Variability")

# Get the density at that depth
rho3d <- map(tsamp,get_3Dslice,ncurl=url,varname='rhopot0',min_zind=1,max_zind=36)

# for each month, extract the mean rho at the depth of max temperature variability
dinds <- match(depth_max_var$depth,dims$z_l$value)
rhomatch <- map2_dbl(rho3d,dinds,\(x,y){
  mean(x[,,y], na.rm=TRUE)
})

tibble(month=1:12,rho=rhomatch) |> 
  ggplot(aes(month,rho))+
  geom_line()+
  scale_x_continuous(breaks=seq_len(12))

# Pick candidate isopycnals
# based on the above, let's do 1024.6-1025.2
target_isopycnals <- 1000+seq(24.6,25.2,by=0.2)

# try with function above
test <- calc_isopyc(url,tidx=1,target_rho=target_isopycnals[1],min_zind = 1,max_zind=36,rast = T)
test2 <- calc_isopyc(url,tidx=1,target_rho=target_isopycnals[2],min_zind = 1,max_zind=36,rast = T)
test3 <- calc_isopyc(url,tidx=1,target_rho=target_isopycnals[3],min_zind = 1,max_zind=36,rast = T)
test4 <- calc_isopyc(url,tidx=1,target_rho=target_isopycnals[4],min_zind = 1,max_zind=36,rast = T)
dev.off()
par(mfrow=c(2,2))
plot(test,main=target_isopycnals[1])
plot(test2,main=target_isopycnals[2])
plot(test3,main=target_isopycnals[3])
plot(test4,main=target_isopycnals[4])

## see how long this would take to run for a year of monthly data and the 4 isopycnals
jobtbl <- crossing(tidx=1:12,target_rho=target_isopycnals) |> 
  mutate(min_zind=1,max_zind=36,rast=TRUE,ncurl=url,jobid=row_number())
rhol <- pmap(jobtbl,calc_isopyc,.progress=T)
# split by isopycnal
rhol_iso <- split(rhol,factor(jobtbl$target_rho)) |> map(rast)
# plot
rhol_iso[[1]] |> panel(main="",col=map.pal("magma",100))

# library(furrr)
# plan(multisession,workers=4)
# rhol <- future_pmap(jobtbl,calc_isopyc,.progress = T)

# Function to use interpolation to track temperature along an isopycnal
# rho <- rho3d[[1]]
# temp <- t3d[[1]]
# isopyc <- 1025.0
# 
# interp_isopyc <- function(rho,temp,isopyc=1025.0){
#   
#   if(!all(dim(rho)==dim(temp))) stop("rho and temp have different dimensions.")
#   
#   lon_i <- dim(rho)[1]
#   lat_j <- dim(rho)[2]
#   ll <- crossing(i=seq_len(lon_i),j=seq_len(lat_j))
#   
#   tint_vec <- pmap(ll,\(i,j){
#     rho_prof <- rho[i, j, ]
#     temp_prof <- temp[i, j, ]
#     # Remove NAs for interpolation
#     valid <- !is.na(rho_prof) & !is.na(temp_prof)
#     
#     if (sum(valid) < 2) return(NA)  # Need at least 2 points
#     else{
#       rho_clean <- rho_prof[valid]
#       temp_clean <- temp_prof[valid]
#       # Check if target density is within range
#       out <- NA
#       if (isopyc >= min(rho_clean) && isopyc <= max(rho_clean)) {
#         # Linear interpolation
#         tint <- approx(x = rho, y = temp, 
#                        xout = isopyc, method = "linear")
#         # return the interpolated temperature value
#         out <- tint$y
#       }}
#       out
#     })
#   tint_mat <- matrix(nrow=dim(temp)[1],ncol=dim(temp)[2],data = tint_vec)
#   tint_mat
# }
#   
# for (i in 1:nlat) {
#   for (j in 1:nlon) {
#     rho_prof <- rho[i, j, ]
#     temp_prof <- temp[i, j, ]
#     
#     # Remove NAs for interpolation
#     valid <- !is.na(rho_prof) & !is.na(temp_prof)
#     
#     if (sum(valid) < 2) next  # Need at least 2 points
#     
#     rho_clean <- rho_prof[valid]
#     temp_clean <- temp_prof[valid]
#     
#     # Check if target density is within range
#     if (target_rho >= min(rho_clean) && target_rho <= max(rho_clean)) {
#       # Linear interpolation (approx is like scipy's interp1d with kind='linear')
#       interp_result <- approx(x = rho_clean, y = temp_clean, 
#                               xout = target_rho, method = "linear")
#       temp_on_isopycnal[i, j] <- interp_result$y
#     }
#   }
# }