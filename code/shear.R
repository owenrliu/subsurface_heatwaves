## Calculate vertical shear
library(terra)
library(here)

dims <- readRDS(here('data','mom6',"dims_meta.rds"))
z <- dims$z_l$value[1:36]            # layer centers (m)
dz <- diff(z)                         # spacing between adjacent centers
zmid <- (z[-1] + z[-length(z)]) / 2   # shear is located midway between layers

# calculate shear between depth layers
shear_month <- function(uo, vo) {     # uo, vo: 36-layer SpatRasters for one month
  k  <- 1:(nlyr(uo) - 1)
  dudz <- (uo[[k + 1]] - uo[[k]]) / dz
  dvdz <- (vo[[k + 1]] - vo[[k]]) / dz
  S <- sqrt(dudz^2 + dvdz^2)          # s^-1
  names(S) <- paste0("z", round(zmid, 1))
  S
}

fu <- here('data','mom6','uo_rotate','uo_rotate_1997.nc')
fv <- here('data','mom6','vo_rotate','vo_rotate_1997.nc')

u <- sds(fu); v <- sds(fv)  # 12 subdatasets = months
S <- lapply(seq_along(u), \(m) shear_month(u[m], v[m]))
S_50_150 <- rast(lapply(S, \(s) mean(s[[zmid >= 50 & zmid < 150]], na.rm = TRUE)))

# all years
walk(1993:2025,\(y){
  fu <- here('data','mom6','uo_rotate',paste0('uo_rotate_',y,'.nc'))
  fv <- here('data','mom6','vo_rotate',paste0('vo_rotate_',y,'.nc'))
  
  u <- sds(fu); v <- sds(fv)  # 12 subdatasets = months
  S <- lapply(seq_along(u), \(m) shear_month(u[m], v[m]))
  S_100_200 <- rast(lapply(S, \(s) mean(s[[zmid >= 100 & zmid < 200]], na.rm = TRUE)))
  writeRaster(S_100_200,filename = here('data','mom6',paste0('shear_100_200_',y,'.tif')))
})

# one big baby
fls <- list.files(here('data','mom6'),full.names = T) |> str_subset("shear_")

shear_all <- map(fls,rast) |> rast()
