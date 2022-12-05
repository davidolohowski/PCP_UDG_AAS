library(tidyverse)
library(spatstat)
library(spatial)
library(sp)
library(sf)
library(rgeos)

# load functions and data
source('src/BDMMH_supp.R')
source('src/BDMMH_col.R')
v_acs <- read_csv('data/v11wfc3_GC.csv')
col <- readRDS('data/v11_wfc3_col.RDS')

# specify study region S
X <- c(1, 4129, 4216, 40)/4400
Y <- c(1, 223, 4382, 4177)/4400
S <- Polygon(cbind(X,Y))
S <- SpatialPolygons(list(Polygons(list(S),'region')))
S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))

pts <- v_acs[,3:4]/4400

plot(pts, asp = 1)

b0_start <- log(30/sf::st_area(sf::st_as_sf(S)))
s0_start <- log(0.2)
UDG_start <- data.frame(N = 4, R = 150/4300, e = 1, n = 1, theta = pi/2, sigma = 0.15)
lc_start <- log(1/sf::st_area(sf::st_as_sf(S)))

# specify prior parameters
UDG_prior <- data.frame(N_mean = log(3), N_sd = 1, R_mean = -3.9, R_sd = 0.5,
                        e_mean = log(1.1), e_sd = 0.3, n_mean = 0, n_sd = 0.75)

b0_prior <- c(b0_start, 0.25)
lc_prior <- c(log(1), 0.2)
s0_prior <- c(0.2, 0.05)

# run model 1
seed <- 632178
set.seed(seed)
system.time(res <- BDM_MH_supperpose(X = pts, S = S, b0_start = b0_start, 
                                     UDG_start = UDG_start, lc_start = lc_start, npix = 90000,
                                     UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, 
                                     M = 100000, norm_gal = FALSE))

# run model 2
set.seed(seed)
system.time(res_col <- BDM_MH_supperpose_color(X = pts, Y = col, S = S, b0_start = b0_start, 
                                               UDG_start = UDG_start, s0_start = s0_start, lc_start = lc_start, npix = 90000,
                                               UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, s0_prior = s0_prior, 
                                               M = 100000, norm_gal = FALSE))

#save results
date <- Sys.Date()
filename <- paste0('data/v11wfc3-pcp-col-',date,'-seed-', seed, '-100k.RDS')
saveRDS(res_col, filename)

date <- Sys.Date()
filename <- paste0('data/v11wfc3-pcp-',date,'-seed-', seed, '-100k.RDS')
saveRDS(res, filename)
