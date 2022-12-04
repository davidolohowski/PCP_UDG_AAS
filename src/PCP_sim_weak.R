library(tidyverse)
library(spatstat)
library(spatial)
library(sp)
library(sf)
library(rgeos)
library(emdbook)
library(swfscMisc)
library(INLA)
library(coda)

source('src/BDMMH_supp.R')
source('src/BDMMH_col.R')
v_acs <- read_csv('data/v11acs_GC.csv')

# specify study region S
X <- c(1, 4095, 4214, 219)/4300
Y <- c(1, 100.5, 4300, 4044)/4300
S <- Polygon(cbind(X,Y))
S <- SpatialPolygons(list(Polygons(list(S),'region')))
S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))


U1 <- c(2628, 1532)/4300
U2 <- c(2517, 3289)/4300

U1_GC <- which(sqrt((v_acs$x/4300 - U1[1])^2 + (v_acs$y/4300 - U1[2])^2) < 0.054)
U2_GC <- which(sqrt((v_acs$x/4300 - U2[1])^2 + (v_acs$y/4300 - U2[2])^2) < 0.018)

mean_U1_col <- mean(v_acs$color[U1_GC])
mean_U2_col <- mean(v_acs$color[U2_GC])
sd_U1_col <- sd(v_acs$color[U1_GC])
sd_U2_col <- sd(v_acs$color[U2_GC])

pts <- v_acs %>%
  select(x, y)/4300

pts <- pts[-c(65, 138),]
plot(pts)

# specify starting values
b0_start <- log(100/sf::st_area(sf::st_as_sf(S)))
s0_start <- log(0.2)
UDG_start <- data.frame(N = 6, R = 150/4300, e = 1, n = 1, theta = pi/2, sigma = 0.15)
lc_start <- log(2/sf::st_area(sf::st_as_sf(S)))
gal_start <- cbind(N = c(log(100), log(250)), R = c(6.56 - log(4300), 6.8 - log(4300)), n = c(0.5, 0.5))

# specify prior parameters
UDG_prior <- data.frame(N_mean = log(6), N_sd = 1, R_mean = -3.9, R_sd = 0.5,
                        e_mean = log(1.1), e_sd = 0.3, n_mean = 0, n_sd = 0.75)

b0_prior <- c(b0_start, 0.25)
lc_prior <- c(log(3), 0.2)
s0_prior <- c(0.2, 0.05)

gal_prior <- data.frame(N_mean = c(log(100), log(250)), N_sd = c(0.25, 0.25),
                        R_mean = c(6.56 - log(4300), 6.8 - log(4300)), R_sd = c(0.25, 0.2),
                        n_mean = c(log(0.5), log(0.5)), n_sd = c(0.5, 0.5))

# fixed parameters for normal galaxies
gal_fix <- data.frame(x = c(1125, 2300)/4300, y = c(1780, 5800)/4300,
                      theta = c(pi/18, pi/6), e = c(1.3166, 1.5))


seed <- 810807
set.seed(seed)
system.time(res <- BDM_MH_supperpose(X = pts, S = S, b0_start = b0_start, gal_fix = gal_fix, gal_start = gal_start, 
                                     UDG_start = UDG_start, lc_start = lc_start, npix = 90000,
                                     gal_prior = gal_prior, UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, 
                                     M = 100000, norm_gal = TRUE))

ncol_dat <- res[[1]][-c(1:10000)]

C_post_r <- lapply(ncol_dat, function(x) x)[seq(1, 90001, by = 5)]
X_nc <- do.call(rbind, C_post_r)
X_nc <- as.data.frame(X_nc)

ggplot(X_nc, aes(x, y)) + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 200, h = c(bw.nrd(X_nc$x), bw.nrd(X_nc$y)))

region <- S@polygons[[1]]@Polygons[[1]]@coords
mesh <- inla.mesh.2d(loc.domain = region, max.edge = c(100/4300, 1000/4300))

spde <- inla.spde2.matern(mesh, alpha = 2)

Amat <- inla.spde.make.A(mesh, loc = as.matrix(pts))
s.index <- inla.spde.make.index(name = 'field', n.spde = spde$n.spde)

seed <- 810807
for (i in 7:100) {
  set.seed(i^3 + i)
  U1_col <- rnorm(6, mean_U1_col, sd_U1_col)
  U2_col <- rnorm(3, mean_U2_col, sd_U2_col)
  
  color <- v_acs$color[-c(65,138)]
  color[c(60, 65, 67, 72, 195, 201)] <- U1_col
  color[c(140, 141, 222)] <- U2_col
  dat_stack <- inla.stack(data = list(color = color),
                          A = list(Amat),
                          effects = list(c(s.index, list(Intercept = 1))), tag = 'data')
  
  Apred <- Amat
  pred_stack <- inla.stack(data = list(color = NA),
                           A = list(Apred),
                           effects = list(c(s.index, list(Intercept = 1))), tag = 'pred')
  
  join_stack <- inla.stack(dat_stack, pred_stack)
  
  form <- color ~ -1 + Intercept + f(field, model = spde)
  
  mod <- inla(form, data = inla.stack.data(join_stack, spde = spde), family = 'gaussian',
              control.predictor = list(A = inla.stack.A(join_stack), compute = T),
              control.compute = list(cpo = T, dic = T))
  
  index.pred <- inla.stack.index(join_stack, "pred")$data
  
  inla_pred <- mod$summary.fitted.values[index.pred, 'mean']
  
  plot(pts$x, color)
  points(pts$x, inla_pred, col = 'red')
  col <- color - inla_pred
  
  # run model 2
  set.seed(seed)
  system.time(res_col <- BDM_MH_supperpose_color(X = pts, Y = col, S = S, b0_start = b0_start, gal_fix = gal_fix, gal_start = gal_start, 
                                                 UDG_start = UDG_start, s0_start = s0_start, lc_start = lc_start, npix = 90000,
                                                 gal_prior = gal_prior, UDG_prior = UDG_prior, b0_prior = b0_prior, lc_prior = lc_prior, s0_prior = s0_prior, 
                                                 M = 100000, norm_gal = TRUE))
  
  col_dat <- res_col[[1]][-c(1:10000)]
  
  C_post_rc <- lapply(col_dat, function(x) x)[seq(1, 90001, by = 5)]
  X_c <- do.call(rbind, C_post_rc)
  X_c <- as.data.frame(X_c)
  
  ggplot(X_c, aes(x, y)) + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 200, h = c(bw.nrd(X_c$x), bw.nrd(X_c$y)))
  
  U1_df <- data.frame(x = 2628, y = 1532)/4300
  U2_df <- data.frame(x = 2517, y = 3289)/4300
  
  U <- rbind(U1_df, U2_df)
  
  ci <- c(0.1, 0.2, 0.3, 0.4)
  dm <- sapply(ci, function(x) detect_metric(U, c(0.018, 0.015), x, X_nc, X_c))
  
  output_file_sim <- paste0('PCP_N=6X3_sim', i, '.RDS')
  saveRDS(dm, output_file_sim)
}
