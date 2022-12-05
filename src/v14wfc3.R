library(tidyverse)
library(spatstat)
library(spatial)
library(sp)
library(sf)
library(rgeos)

# load functions and data
source('src/BDMMH_supp.R')
source('src/BDMMH_col.R')
v_acs <- read_csv('data/v14wfc3_GC.csv')
col <- readRDS('data/v14_wfc3_col.RDS')

# specify study region S
X <- c(1, 4129, 4216, 40)/4400
Y <- c(1, 223, 4382, 4177)/4400
S <- Polygon(cbind(X,Y))
S <- SpatialPolygons(list(Polygons(list(S),'region')))
S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))

# extract point patterns and scale it to [0,1]x[0,1]
pts <- v_acs[,3:4]/4400

# specify starting values
b0_start <- log(30/sf::st_area(sf::st_as_sf(S)))
s0_start <- log(0.2)
UDG_start <- data.frame(N = 5, R = 150/4300, e = 1, n = 1, theta = pi/2, sigma = 0.15)
lc_start <- log(1/sf::st_area(sf::st_as_sf(S)))
gal_start <- cbind(N = c(log(10)), R = c(4.9 - log(4400)), n = c(0.5))

# specify prior parameters
UDG_prior <- data.frame(N_mean = log(3), N_sd = 1, R_mean = -3.9, R_sd = 0.5,
                        e_mean = log(1.1), e_sd = 0.3, n_mean = 0, n_sd = 0.75)

b0_prior <- c(b0_start, 0.25)
lc_prior <- c(log(1), 0.2)
s0_prior <- c(0.2, 0.05)

gal_prior <- data.frame(N_mean = c(log(10)), N_sd = c(0.25),
                        R_mean = c(4.9 - log(4400)), R_sd = c(0.25),
                        n_mean = c(log(0.5)), n_sd = c(0.5))

# fixed parameters for normal galaxies
gal_fix <- data.frame(x = c(1977)/4400, y = c(2962)/4400,
                      theta = c(pi/3), e = c(1.44))


# run model 1
seed <- 209886
set.seed(seed)
system.time(res <- BDM_MH_supperpose(X = pts, S = S, b0_start = b0_start, gal_fix = gal_fix, gal_start = gal_start, 
                                     UDG_start = UDG_start, lc_start = lc_start, npix = 90000,
                                     gal_prior = gal_prior, UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, 
                                     M = 100000, norm_gal = TRUE))

# run model 2
set.seed(seed)
system.time(res_col <- BDM_MH_supperpose_color(X = pts, Y = col, S = S, b0_start = b0_start, gal_fix = gal_fix, gal_start = gal_start, 
                                               UDG_start = UDG_start, s0_start = s0_start, lc_start = lc_start, npix = 90000,
                                               gal_prior = gal_prior, UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, s0_prior = s0_prior, 
                                               M = 100000, norm_gal = TRUE))

#save results
date <- Sys.Date()
filename <- paste0('data/v14wfc3-pcp-col-',date,'-seed-', seed, '-100k.RDS')
saveRDS(res_col, filename)

date <- Sys.Date()
filename <- paste0('data/v14wfc3-pcp-',date,'-seed-', seed, '-100k.RDS')
saveRDS(res, filename)


