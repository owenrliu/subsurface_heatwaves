# Process CORA (interpolated observational) data
# Downloaded on 10/01/2026
library(here)
library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(tictoc)
library(colorspace)

# raster
datr <- rast(here('data','CORA','cora_temperature.nc'),subds="TEMP")

# we want to split this by time layers, but let's make sure they're in the right order
head(names(datr));head(time(datr));head(depth(datr))
# depth varies fastest
nd <- length(unique(depth(datr))) # number of depth levels
nt <- nlyr(datr) / nd             # number of time steps

# Reorganize into a list by time
tind_df <- tibble(lyr=seq_len(nlyr(datr))) |> 
  mutate(dind=rep(1:nd,nt),
         depth=depth(datr),
         tind=rep(1:nt,each=nd),
         time=time(datr)) |> 
  mutate(year=year(time),month=month(time))

by_time <- map(seq_len(nt), \(t) {
  inds <- tind_df |> filter(tind==t) |> pull(lyr)
  datr[[inds]]
})
monyr <- paste(month.name[month(time(datr))],year(time(datr)),sep = "_") |> unique()

# squish back into a SpatRasterDataset
dats <- sds(by_time)            # SpatRasterDataset: s[d] is a 3D block for time t
names(dats) <- monyr

dvals <- tibble(value=sort(unique(depth(datr)))) |> 
  mutate(dbot=value+(lead(value)-value)/2) |> 
  mutate(dbot=coalesce(dbot,1000)) |>
  # layer thickness
  mutate(dz=diff(c(0,dbot)))

# Slice and calculate OHC
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

test <- dats[[78]]
ohcts <- app(test,make_ohc,minz=0,maxz=50)
plot(ohcts)

# apply function to all layers
cora_ohc_z_0_50 <- imap(dats,\(x,n){
  tic(paste("Layer",n,"processed."))
  out <- app(x,make_ohc,minz=0,maxz=50)
  toc()
  out
}) |> rast()

cora_ohc_z_50_100 <- imap(dats,\(x,n){
  tic(paste("Layer",n,"processed."))
  out <- app(x,make_ohc,minz=50,maxz=100)
  toc()
  out
}) |> rast()

cora_ohc_z_100_250 <- imap(dats,\(x,n){
  tic(paste("Layer",n,"processed."))
  out <- app(x,make_ohc,minz=100,maxz=250)
  toc()
  out
}) |> rast()

cora_ohc_z_250_500 <- imap(dats,\(x,n){
  tic(paste("Layer",n,"processed."))
  out <- app(x,make_ohc,minz=250,maxz=500)
  toc()
  out
}) |> rast()

cora_ohc_z_500_1000 <- imap(dats,\(x,n){
  tic(paste("Layer",n,"processed."))
  out <- app(x,make_ohc,minz=500,maxz=1000)
  toc()
  out
}) |> rast()

time(cora_ohc_z_0_50) <- time(cora_ohc_z_50_100) <- time(cora_ohc_z_100_250) <- time(cora_ohc_z_250_500) <- time(cora_ohc_z_500_1000) <- unique(tind_df$time)

writeRaster(cora_ohc_z_0_50,here('output','cora_ohc_0_50_ts.tif'),overwrite=T)
writeRaster(cora_ohc_z_50_100,here('output','cora_ohc_50_100_ts.tif'),overwrite=T)
writeRaster(cora_ohc_z_100_250,here('output','cora_ohc_100_250_ts.tif'),overwrite=T)
writeRaster(cora_ohc_z_250_500,here('output','cora_ohc_250_500_ts.tif'),overwrite=T)
writeRaster(cora_ohc_z_500_1000,here('output','cora_ohc_500_1000_ts.tif'),overwrite=T)
