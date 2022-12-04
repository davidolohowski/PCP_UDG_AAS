library(tidyverse)
library(spatstat)
library(spatial)
library(sp)
library(sf)
library(rgeos)
library(plotly)

# load functions and data
source('src/BDMMH_supp.R')
source('src/BDMMH_col.R')
v_acs <- read_csv('data/v12acs_GC.csv')
col <- readRDS('data/v12_acs_col.RDS')

# specify study region S
X <- c(1, 4095, 4214, 219)/4300
Y <- c(1, 100.5, 4300, 4044)/4300
S <- Polygon(cbind(X,Y))
S <- SpatialPolygons(list(Polygons(list(S),'region')))
S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))

# extract point patterns and scale it to [0,1]x[0,1]
pts <- v_acs[,3:4]/4300

plot(pts, asp = 1)

# specify starting values
b0_start <- log(50/sf::st_area(sf::st_as_sf(S)))
s0_start <- log(0.2)
UDG_start <- data.frame(N = 3, R = 150/4300, e = 1, n = 1, theta = pi/2, sigma = 0.15)
lc_start <- log(1/sf::st_area(sf::st_as_sf(S)))
gal_start <- cbind(N = c(log(10), log(40)), R = c(4.4 - log(4300), 5.5 - log(4300)), n = c(0.5, 0.5))

# specify prior parameters
UDG_prior <- data.frame(N_mean = log(3), N_sd = 1, R_mean = -3.9, R_sd = 0.5,
                        e_mean = log(1.1), e_sd = 0.3, n_mean = 0, n_sd = 0.75)

b0_prior <- c(b0_start, 0.25)
lc_prior <- c(log(1), 0.2)
s0_prior <- c(0.2, 0.05)

gal_prior <- data.frame(N_mean = c(log(10), log(40)), N_sd = c(0.25, 0.25),
                        R_mean = c(4.4 - log(4300), 5.5 - log(4300)), R_sd = c(0.25, 0.25),
                        n_mean = c(log(0.5), log(0.5)), n_sd = c(0.5, 0.5))

# fixed parameters for normal galaxies
gal_fix <- data.frame(x = c(3108, 2514)/4300, y = c(3377, -10)/4300,
                      theta = c(5*pi/180, pi/6), e = c(1.14, 2))


# run model 1
seed <- 513268
set.seed(seed)
system.time(res <- BDM_MH_supperpose(X = pts, S = S, b0_start = b0_start, gal_fix = gal_fix, gal_start = gal_start, 
                                     UDG_start = UDG_start, lc_start = lc_start, npix = 90000,
                                     gal_prior = gal_prior, UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, 
                                     M = 100000, norm_gal = T))

set.seed(seed)
system.time(res_col <- BDM_MH_supperpose_color(X = pts, Y = col, S = S, b0_start = b0_start, gal_fix = gal_fix, gal_start = gal_start, 
                                               UDG_start = UDG_start, s0_start = s0_start, lc_start = lc_start, npix = 90000,
                                               gal_prior = gal_prior, UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, s0_prior = s0_prior, 
                                               M = 100000, norm_gal = T))

#save results
date <- Sys.Date()
filename <- paste0('data/v12acs-pcp-col-',date,'-seed-', seed, '-100k.RDS')
saveRDS(res_col, filename)

date <- Sys.Date()
filename <- paste0('data/v12acs-pcp-',date,'-seed-', seed, '-100k.RDS')
saveRDS(res, filename)
