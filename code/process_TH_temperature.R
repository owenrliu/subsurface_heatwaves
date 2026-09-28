# process trinidad head temperature
library(rerddap)
library(here)
library(tidyverse)

# Query the ERDDAP dataset
dat <- tabledap(
  "swfscTrinidadCTD",
  fields = c("station", "time", "latitude", "longitude","depth", "temperature"),
  url = "https://oceanview.pfeg.noaa.gov/erddap/")

dat <- dat |> 
  mutate(stationfct=factor(station)) |> 
  mutate(year=year(time),month=month(time),date=date(time)) |> 
  mutate(depthbin=cut(depth,breaks=seq(0,550,by=10)))

dat |> 
  ggplot(aes(depth,temperature))+
  geom_point()
