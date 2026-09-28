# Process California Underwater Glider Network data
# downloaded from
# https://spraydata.ucsd.edu/projects/cugn/#about-project-dataset

library(here)
library(tidyverse)
library(ncdf4)
library(terra)
library(viridis)
library(colorspace)

# nc <- nc_open(here('data','CUGN_line_90.nc'))
nc <- nc_open(here('data','anomaly_z_90.nc'))
# variable names as a table of definitions
vardefs <- map(names(nc$var),\(x){
  ncatt_get(nc,x) |> 
    as_tibble() |> 
    mutate(varname=x)
}) |> bind_rows()

# dimensions as a list of indexed lookup tables
dims <- imap(nc$dim, \(x,n){
  tibble(value=c(x$vals),name=n,units=x$units) |> 
    mutate(index=row_number())
}) |> set_names(names(nc$dim))

# all the data (time x depth x distance)
dat <- ncvar_get(nc,varid="temperature")

# data for 100km offshore
dat100 <- dat[,,21]

#lat/lon reference for dive midpoint
# lats <- ncvar_get(nc,varid="lat") |> c()
# lons <- ncvar_get(nc,varid="lon") |> c()

# depths
depthref <- dims$depth

# times
times <- dims$time |> 
  mutate(date=as_date(value,origin="1970-01-01"))

# make a table
df <- crossing(depth=depthref$value,date=times$date) |> 
  mutate(temp=c(dat100)) |> 
  mutate(ym=round_date(date,"month"))

# measure distance from coast
# coast <- ne_states(country='United States of America',returnclass = 'sf') |> 
#   filter(name %in% c('California'))
# 
# ll <- df |> select(lon,lat) |> st_as_sf(coords=c('lon','lat'),crs=4326)
# dists <- st_distance(ll,coast)
# df$distcoast_m <- dists[,1] |> as.numeric()
# df$distcoast_km <- df$distcoast_m/1000

test <- df |> 
  group_by(ym,depth) |> 
  summarise(meananom=mean(temp,na.rm=T),.groups='drop')

# NAs?
testNA <- test |> group_by(depth) |> summarise(isNA=sum(is.na(meananom)))
testNAdate <- test |> group_by(ym) |> summarise(isNA=sum(is.na(meananom)))
# big break in data availability after 280m depth

test |> 
  filter(depth <= 280) |> 
  ggplot(aes(ym,meananom,color=as.factor(depth)))+geom_point()

p_

p_anom_hov_100k <- test |> 
  filter(depth <= 280,ym>"2006-10-01",ym<"2026-09-01") |>
  ggplot(aes(ym,depth,fill=meananom,col=meananom))+
  geom_tile()+
  scale_y_reverse(expand=c(0,0),breaks=seq(50,250,by=50))+
  scale_x_date(expand=c(0,0),date_breaks = '2 years',date_labels = "%Y")+
  scale_fill_continuous_diverging(palette = "Blue-Red 2", mid = 0)+
  scale_color_continuous_diverging(palette = "Blue-Red 2", mid = 0)+
  labs(x="Year",y="Depth (m)",fill="Anomaly",color="Anomaly",
       title="CUGN Line 90 Temperature Anomaly\nat 100km offshore")

p_anom_hov_100k
ggsave(here('plots','cugn_line90_100km_hov.png'),p_anom_hov_100k,w=7)

# mean across depth bins
p_ts_depth <- test |> 
  mutate(depthcat=case_when(
    depth<=50 ~ "0-50",
    depth>50&depth<=100 ~ "50-100",
    depth>100&depth<=250 ~ "100-250"
  )) |> 
  mutate(depthcat=factor(depthcat,levels=c("0-50","50-100","100-250"))) |> 
  group_by(ym,depthcat) |> 
  summarise(meananom=mean(meananom,na.rm=T),.groups='drop') |> 
  filter(!is.na(depthcat)) |> 
  ggplot(aes(ym,meananom,color=depthcat,linetype=depthcat))+
  geom_hline(yintercept=0,linetype=2)+
  geom_line(linewidth=1)+
  scale_color_manual(values=viridis_pal(option="G",end=0.7)(3))+
  labs(x="Year",y="Mean Anomaly",color='Depth',linetype="Depth")
p_ts_depth

ggsave(here('plots','cugn_line90_100km_ts.png'),p_ts_depth,w=10)
