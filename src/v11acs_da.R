library(tidyverse)
library(modeest)
library(coda)
library(stableGR)
library(RColorBrewer)
library(reshape2)
library(emdbook)
library(swfscMisc)
library(tikzDevice)
library(raster)
library(spatial)
library(expint)
library(sf)
library(gplots)

options(warn = -1)
source('src/help_functions.R')
v_acs <- read_csv('data/v11acs_GC.csv')
col <- readRDS('data/v11_acs_col.RDS')
pts <- v_acs[,3:4]/4300

# specify study region S
X <- c(1, 4095, 4214, 219)/4300
Y <- c(1, 100.5, 4300, 4044)/4300
S <- Polygon(cbind(X,Y))
S <- SpatialPolygons(list(Polygons(list(S),'region')))
S <- SpatialPolygonsDataFrame(S, data.frame(id = S@polygons[[1]]@ID, row.names = S@polygons[[1]]@ID))

gal_fix <- data.frame(x = c(1125, 2300)/4300, y = c(1780, 5800)/4300,
                      theta = c(pi/18, pi/6), e = c(1.3166, 1.5))

# chain for model 1
res <- readRDS('data/v11acs_res/v11acs-pcp-2022-07-16-seed-810807-100k.RDS')

# chain for model 2
res_col <- readRDS('data/v11acs_res/v11acs-pcp-col-2022-07-16-seed-810807-100k.RDS')

# extract parameters to compute P(Z != 0 | X) (first 10000 as bur-in)
sub_res <- res[[1]]
sub_fix <- res[[2]]

b0 <- exp(sub_fix[,'b0'])
gal1_par <- exp(sub_fix[,3:5])
gal2_par <- exp(sub_fix[,6:8])

colnames(gal1_par) <- c('N', 'R', 'n')
colnames(gal2_par) <- c('N', 'R', 'n')

gal1_par <- as.data.frame(gal1_par)
gal2_par <- as.data.frame(gal2_par)

gal1_par$x <- gal_fix$x[1]
gal1_par$y <- gal_fix$y[1]
gal1_par$theta <- gal_fix$theta[1]
gal1_par$e <- gal_fix$e[1]

gal2_par$x <- gal_fix$x[2]
gal2_par$y <- gal_fix$y[2]
gal2_par$theta <- gal_fix$theta[2]
gal2_par$e <- gal_fix$e[2]

gal_par <- list()

for(j in 1:100001){
  gal_par[[j]] <- rbind(gal1_par[j,], gal2_par[j,])
}

prob <- matrix(0, 262, 100001)
for (j in 1:100001) {
  prob[,j] <- PZ_X(pts, sub_res[[j]], gal_par[[j]], b0[j]) 
}
prob_ncol <- t(prob)

# extract parameters to compute P(Z |!= 0 | X,V) (first 10000 as bur-in)
sub_res_col <- res_col[[1]]
sub_fix_col <- res_col[[2]]

b0_col <- exp(sub_fix_col[,'b0'])
s0_col <- exp(sub_fix_col[,'sigma0'])
gal1_par_col <- exp(sub_fix_col[,4:6])
gal2_par_col <- exp(sub_fix_col[,7:9])

colnames(gal1_par_col) <- c('N', 'R', 'n')
colnames(gal2_par_col) <- c('N', 'R', 'n')

gal1_par_col <- as.data.frame(gal1_par_col)
gal2_par_col <- as.data.frame(gal2_par_col)

gal1_par_col$x <- gal_fix$x[1]
gal1_par_col$y <- gal_fix$y[1]
gal1_par_col$theta <- gal_fix$theta[1]
gal1_par_col$e <- gal_fix$e[1]

gal2_par_col$x <- gal_fix$x[2]
gal2_par_col$y <- gal_fix$y[2]
gal2_par_col$theta <- gal_fix$theta[2]
gal2_par_col$e <- gal_fix$e[2]

gal_par <- list()

for(j in 1:100001){
  gal_par[[j]] <- rbind(gal1_par_col[j,], gal2_par_col[j,])
}

prob <- matrix(0, 262, 100001)
for (j in 1:100001) {
  prob[,j] <- PZ_XY(pts, col, sub_res_col[[j]], gal_par[[j]], b0_col[j], s0_col[j]) 
}

prob_col <- t(prob)

prob_ncol <- prob_ncol[-c(1:10001),]
prob_col <- prob_col[-c(1:10001),]

# R hat statsitics for Model 1 and Model 2
PSRF_ncol <- stable.GR(prob_ncol, multivariate = F)
PSRF_col <- stable.GR(prob_col, multivariate = F)

# posterior density plots
ncol_mod <- as.data.frame(exp(res[[2]][-c(1:10000),]))
prior_df_ncol <- data.frame(value = c(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 
                                      seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100),
                                      seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100),
                                      seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100),
                                      seq(min(ncol_mod$n2), max(ncol_mod$n2), length= 100),
                                      seq(min(ncol_mod$N2), max(ncol_mod$N2), length= 100),
                                      seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100),
                                      seq(min(ncol_mod$R2), max(ncol_mod$R2), length= 100)),
                            density = c(dlnorm(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 4.711201, 0.25),
                                        dlnorm(seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100), log(3), 0.2),
                                        dlnorm(seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100), log(0.5), 0.5),
                                        dlnorm(seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100), log(100), 0.25),
                                        dlnorm(seq(min(ncol_mod$n2), max(ncol_mod$n2), length= 100), log(0.5), 0.5),
                                        dlnorm(seq(min(ncol_mod$N2), max(ncol_mod$N2), length= 100), log(250), 0.25),
                                        dlnorm(seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100), 6.56 - log(4300), 0.25),
                                        dlnorm(seq(min(ncol_mod$R2), max(ncol_mod$R2), length= 100), 6.8 - log(4300), 0.2)),
                            key = c(rep('b0', 100), rep('lc', 100), rep('n1', 100), rep('N1', 100), rep('n2', 100), rep('N2', 100),
                                    rep('R1', 100), rep('R2', 100)))

col_mod <- as.data.frame(exp(res_col[[2]][-c(1:10000),]))
prior_df_col <- data.frame(value = c(seq(min(col_mod$b0), max(col_mod$b0), length= 100), 
                                     seq(min(col_mod$lc), max(col_mod$lc), length= 100),
                                     seq(min(col_mod$n1), max(col_mod$n1), length= 100),
                                     seq(min(col_mod$N1), max(col_mod$N1), length= 100),
                                     seq(min(col_mod$n2), max(col_mod$n2), length= 100),
                                     seq(min(col_mod$N2), max(col_mod$N2), length= 100),
                                     seq(min(col_mod$R1), max(col_mod$R1), length= 100),
                                     seq(min(col_mod$R2), max(col_mod$R2), length= 100),
                                     seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100)),
                           density = c(dlnorm(seq(min(col_mod$b0), max(col_mod$b0), length= 100), 4.711201, 0.25),
                                       dlnorm(seq(min(col_mod$lc), max(col_mod$lc), length= 100), log(3), 0.2),
                                       dlnorm(seq(min(col_mod$n1), max(col_mod$n1), length= 100), log(0.5), 0.5),
                                       dlnorm(seq(min(col_mod$N1), max(col_mod$N1), length= 100), log(100), 0.25),
                                       dlnorm(seq(min(col_mod$n2), max(col_mod$n2), length= 100), log(0.5), 0.5),
                                       dlnorm(seq(min(col_mod$N2), max(col_mod$N2), length= 100), log(250), 0.25),
                                       dlnorm(seq(min(col_mod$R1), max(col_mod$R1), length= 100), 6.56 - log(4300), 0.25),
                                       dlnorm(seq(min(col_mod$R2), max(col_mod$R2), length= 100), 6.8 - log(4300), 0.2),
                                       dgamma(seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100), 0.2, 0.05)),
                           key = c(rep('b0', 100), rep('lc', 100), rep('n1', 100), rep('N1', 100), rep('n2', 100), rep('N2', 100),
                                   rep('R1', 100), rep('R2', 100), rep('sigma0', 100)))

ncol.labs <- c('$\\beta_0$', '$\\lambda_c$', '$n^1_g$', '$N_g^1$', '$n_2^g$', '$N^2_g$', '$R_g^1$', '$R_g^2$')
names(ncol.labs) <- c('b0', 'lc', 'n1', 'N1', 'n2', 'N2', 'R1', 'R2')

col.labs <- c('$\\beta_0$', '$\\lambda_c$', '$n^1_g$', '$N_g^1$', '$n_2^g$', '$N^2_g$', '$R_g^1$', '$R_g^2$', '$\\sigma_0$')
names(col.labs) <- c('b0', 'lc', 'n1', 'N1', 'n2', 'N2', 'R1', 'R2', 'sigma0')

tikz(file = "v11acs_post_ncol.tex", standAlone=T, width = 9, height = 9)
ggplot(gather(ncol_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_ncol, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol.labs))
dev.off()
system('pdflatex v11acs_post_ncol.tex')

tikz(file = "v11acs_post_col.tex", standAlone=T, width = 9, height = 9)
ggplot(gather(col_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_col, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col.labs))
dev.off()
system('pdflatex v11acs_post_col.tex')

# whether there are UDGs? (Yes)
ncol_dat <- res[[1]][-c(1:10000)]
col_dat <- res_col[[1]][-c(1:10000)]

C_post_r <- lapply(ncol_dat, function(x) x)
X_nc <- do.call(rbind, C_post_r)
X_nc <- as.data.frame(X_nc)

C_post_rc <- lapply(col_dat, function(x) x)
X_c <- do.call(rbind, C_post_rc)
X_c <- as.data.frame(X_c)

N_nc <- 1*(unlist(lapply(ncol_dat, nrow)) == 0)
N_c <- 1*(unlist(lapply(col_dat, nrow)) == 0)

sum(N_nc)/90001
sum(N_c)/90001

p0_nc <- numeric(100)
p0_c <- numeric(100)

pu_nc <- numeric(100)
pu_c <- numeric(100)

hist_res_nc <- hist2d(X_nc$x, X_nc$y, nbins = c(7, 7))
hist_res_c <- hist2d(X_c$x, X_c$y, nbins = c(7, 7))

max(hist_res_nc$counts)/nrow(X_nc)*(1 - sum(N_nc)/90001)
max(hist_res_c$counts)/nrow(X_c)*(1 - sum(N_c)/90001)

set.seed(113409)
for (i in 1:100) {
  block_idx <- sample(1:720, 720, replace = T)
  rs_block <- as.vector(sapply(block_idx, function(x) c((1+125*(x-1)):(x*125))))
  X_nc_j <- ncol_dat[rs_block]
  X_nc_j <- do.call(rbind, X_nc_j)
  X_nc_j <- as.data.frame(X_nc_j)
  X_c_j <- col_dat[rs_block]
  X_c_j <- do.call(rbind, X_c_j)
  X_c_j <- as.data.frame(X_c_j)
  N_nc_j <- 1*(unlist(lapply(ncol_dat[rs_block], nrow)) == 0)
  N_c_j <- 1*(unlist(lapply(col_dat[rs_block], nrow)) == 0)
  p0_nc[i] <- sum(N_nc_j)/90001
  p0_c[i] <- sum(N_c_j)/90001
  hist_res_nc_j <- hist2d(X_nc_j$x, X_nc_j$y, nbins = c(7,7))
  hist_res_c_j <- hist2d(X_c_j$x, X_c_j$y, nbins = c(7,7))
  pu_nc[i] <- max(hist_res_nc_j$counts)/nrow(X_nc_j)*(1 - sum(N_nc_j)/90001)
  pu_c[i] <- max(hist_res_c_j$counts)/nrow(X_c_j)*(1 - sum(N_c_j)/90001)
}

# Posterior distributions regarding UDGs
U1 <- data.frame(x = 2628, y = 1532)/4300
U2 <- data.frame(x = 2517, y = 3289)/4300

U <- rbind(U1, U2)

df_U1_ncol <- do.call(rbind, res[[1]][-c(1:10001)]) %>%
  filter((x - U1$x)^2 + (y - U1$y)^2 <= 0.01^2) %>%
  dplyr::select(N, R, n)

df_U2_ncol <- do.call(rbind, res[[1]][-c(1:10001)]) %>%
  filter((x - U2$x)^2 + (y - U2$y)^2 <= 0.01^2)%>%
  dplyr::select(N, R, n)

df_U1_col <- do.call(rbind, res_col[[1]][-c(1:10001)]) %>%
  filter((x - U1$x)^2 + (y - U1$y)^2 <= 0.01^2)%>%
  dplyr::select(N, R, n, sigma)

df_U2_col <- do.call(rbind, res_col[[1]][-c(1:10001)]) %>%
  filter((x - U2$x)^2 + (y - U2$y)^2 <= 0.01^2)%>%
  dplyr::select(N, R, n, sigma)

prior_df_ncol_U1 <- data.frame(value = c(seq(min(df_U1_ncol$N), max(df_U1_ncol$N), length= 100), 
                                      seq(min(df_U1_ncol$R), max(df_U1_ncol$R), length= 100),
                                      seq(min(df_U1_ncol$n), max(df_U1_ncol$n), length= 100)),
                            density = c(dlnorm(seq(min(df_U1_ncol$N), max(df_U1_ncol$N), length= 100), log(6), 1),
                                        dlnorm(seq(min(df_U1_ncol$R), max(df_U1_ncol$R), length= 100), -3.9, 0.5),
                                        dlnorm(seq(min(df_U1_ncol$n), max(df_U1_ncol$n), length= 100), 0, 0.75)),
                            key = c(rep('N', 100), rep('R', 100), rep('n', 100)))

ncol_U1_labs <- c('$N_u^1$', '$R_u^1$', '$n_u^1$')
names(ncol_U1_labs) <- c('N', 'R', 'n')

tikz(file = "v11acs_U1_post_ncol.tex", standAlone=T, width = 9, height = 3)
ggplot(gather(df_U1_ncol), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_ncol_U1, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol_U1_labs))
dev.off()
system('pdflatex v11acs_U1_post_ncol.tex')

prior_df_ncol_U2 <- data.frame(value = c(seq(min(df_U2_ncol$N), max(df_U2_ncol$N), length= 100), 
                                         seq(min(df_U2_ncol$R), max(df_U2_ncol$R), length= 100),
                                         seq(min(df_U2_ncol$n), max(df_U2_ncol$n), length= 100)),
                               density = c(dlnorm(seq(min(df_U2_ncol$N), max(df_U2_ncol$N), length= 100), log(6), 1),
                                           dlnorm(seq(min(df_U2_ncol$R), max(df_U2_ncol$R), length= 100), -3.9, 0.5),
                                           dlnorm(seq(min(df_U2_ncol$n), max(df_U2_ncol$n), length= 100), 0, 0.75)),
                               key = c(rep('N', 100), rep('R', 100), rep('n', 100)))

ncol_U2_labs <- c('$N_u^2$', '$R_u^2$', '$n_u^2$')
names(ncol_U2_labs) <- c('N', 'R', 'n')

tikz(file = "v11acs_U2_post_ncol.tex", standAlone=T, width = 9, height = 3)
ggplot(gather(df_U2_ncol), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) +  
  geom_line(data = prior_df_ncol_U2, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol_U2_labs))
dev.off()
system('pdflatex v11acs_U2_post_ncol.tex')

prior_df_col_U1 <- data.frame(value = c(seq(min(df_U1_col$N), max(df_U1_col$N), length= 100), 
                                         seq(min(df_U1_col$R), max(df_U1_col$R), length= 100),
                                         seq(min(df_U1_col$n), max(df_U1_col$n), length= 100),
                                         seq(min(df_U1_col$sigma), max(df_U1_col$sigma), length= 100)),
                               density = c(dlnorm(seq(min(df_U1_col$N), max(df_U1_col$N), length= 100), log(6), 1),
                                           dlnorm(seq(min(df_U1_col$R), max(df_U1_col$R), length= 100), -3.9, 0.5),
                                           dlnorm(seq(min(df_U1_col$n), max(df_U1_col$n), length= 100), 0, 0.75),
                                           5/gamma(0.05)*gammainc(0.05 - 1, seq(min(df_U1_col$sigma), max(df_U1_col$sigma), length= 100))),
                               key = c(rep('N', 100), rep('R', 100), rep('n', 100), rep('sigma', 100)))

col_U1_labs <- c('$N_u^1$', '$R_u^1$', '$n_u^1$', '$\\sigma_1$')
names(col_U1_labs) <- c('N', 'R', 'n', 'sigma')

tikz(file = "v11acs_U1_post_col.tex", standAlone=T, width = 6, height = 6)
ggplot(gather(df_U1_col), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_col_U1, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col_U1_labs))
dev.off()
system('pdflatex v11acs_U1_post_col.tex')

prior_df_col_U2 <- data.frame(value = c(seq(min(df_U2_col$N), max(df_U2_col$N), length= 100), 
                                         seq(min(df_U2_col$R), max(df_U2_col$R), length= 100),
                                         seq(min(df_U2_col$n), max(df_U2_col$n), length= 100),
                                        seq(min(df_U2_col$sigma), max(df_U2_col$sigma), length= 100)),
                               density = c(dlnorm(seq(min(df_U2_col$N), max(df_U2_col$N), length= 100), log(6), 1),
                                           dlnorm(seq(min(df_U2_col$R), max(df_U2_col$R), length= 100), -3.9, 0.5),
                                           dlnorm(seq(min(df_U2_col$n), max(df_U2_col$n), length= 100), 0, 0.75),
                                           5/gamma(0.05)*gammainc(0.05 - 1, seq(min(df_U1_col$sigma), max(df_U1_col$sigma), length= 100))),
                               key = c(rep('N', 100), rep('R', 100), rep('n', 100), rep('sigma', 100)))

col_U2_labs <- c('$N_u^2$', '$R_u^2$', '$n_u^2$', '$\\sigma_2$')
names(col_U2_labs) <- c('N', 'R', 'n', 'sigma')

tikz(file = "v11acs_U2_post_col.tex", standAlone=T, width = 6, height = 6)
ggplot(gather(df_U2_col), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_col_U2, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col_U2_labs))
dev.off()
system('pdflatex v11acs_U2_post_col.tex')


# density and plots
ncol_dat <- res[[1]][-c(1:10001)][seq(1, 90000, 5)]
col_dat <- res_col[[1]][-c(1:10001)][seq(1, 90000, 5)]

C_post_r <- lapply(ncol_dat, function(x) x)
X_nc <- do.call(rbind, C_post_r)
X_nc <- as.data.frame(X_nc)

C_post_rc <- lapply(col_dat, function(x) x)
X_c <- do.call(rbind, C_post_rc)
X_c <- as.data.frame(X_c)

# Performance metric comparison
ci <- c(0.1, 0.2, 0.3, 0.4)
dm <- sapply(ci, function(x) detect_metric(U, c(0.018, 0.015), x, X_nc, X_c))
acc_bs <- vector('list', length = 4)
prec_bs <- vector('list', length = 4)
FP_bs <- vector('list', length = 4)
set.seed(277098)
for(i in 1:4){
  acc_bs[[i]] <- numeric(100)
  prec_bs[[i]] <- numeric(100)
  FP_bs[[i]] <- numeric(100)
  for (j in 1:100) {
    block_idx <- sample(1:720, 720, replace = T)
    rs_block <- as.vector(sapply(block_idx, function(x) c((1+25*(x-1)):(x*25))))
    X_nc_j <- C_post_r[rs_block]
    X_nc_j <- do.call(rbind, X_nc_j)
    X_nc_j <- as.data.frame(X_nc_j)
    X_c_j <- C_post_rc[rs_block]
    X_c_j <- do.call(rbind, X_c_j)
    X_c_j <- as.data.frame(X_c_j)
    pm <- detect_metric(U, c(0.018, 0.015), ci[i], X_nc_j, X_c_j)
    acc_bs[[i]][j] <- pm[[1]]
    FP_bs[[i]][j] <- pm[[2]]
    prec_bs[[i]][j] <- pm[[3]]
  }
}

saveRDS(list(dm, acc_bs, FP_bs, prec_bs), 'data/v11acs_res/v11acs_pm.RDS')

# Posterior density plots for X_c
M1_kde <- MASS::kde2d(X_nc$x, X_nc$y, 
                      h = c(bw.nrd(X_nc$x), bw.nrd(X_nc$y)), 
                      n = 1000)

M1_dens <- bind_cols(expand.grid(M1_kde$x, M1_kde$y), z = as.vector(M1_kde$z))
colnames(M1_dens) <- c('X','Y','z')

M2_kde <- MASS::kde2d(X_c$x, X_c$y, 
                      h = c(bw.nrd(X_c$x), bw.nrd(X_c$y)), 
                      n = 1000)

M2_dens <- bind_cols(expand.grid(M2_kde$x, M2_kde$y), z = as.vector(M2_kde$z))
colnames(M2_dens) <- c('X','Y','z')

dens <- bind_rows(M1_dens, M2_dens) %>%
  mutate(Model = c(rep('Model 1', nrow(M1_dens)), rep('Model 2', nrow(M2_dens))))

HPD_dens <- rbind(dm[4,][[1]], dm[4,][[2]], dm[4,][[3]], dm[4,][[4]],
                  dm[5,][[1]], dm[5,][[2]], dm[5,][[3]], dm[5,][[4]])

tikz(file = "v11acs_res.tex", standAlone=T, width = 6, height = 3)
ggplot() + geom_raster(data = dens, aes(X, Y, fill = z)) + 
  geom_sf(data = HPD_dens, aes(color = prob), fill = NA, size = 0.1) +
  coord_sf() +  theme_minimal() + xlim(c(0,1)) + ylim(c(0,1)) +
  facet_grid(~Model) + scale_fill_viridis_c(option = 'D', name = 'Density') + 
  scale_color_viridis_c(option = 'A', direction = -1, name = '$1-\\alpha$') +
  xlab('X') + ylab('Y')
dev.off()
system('pdflatex v11acs_res.tex')

WUDG88_M1_dens <- M1_dens %>%
  filter(X > as.numeric(U1[1]) - 0.018*4 & X < as.numeric(U1[1]) + 0.018*4 &
           Y > as.numeric(U1[2]) - 0.018*4 & Y < as.numeric(U1[2]) + 0.018*4)

WUDG89_M1_dens <- M1_dens %>%
  filter(X > as.numeric(U2[1]) - 0.015*4 & X < as.numeric(U2[1]) + 0.015*4 &
           Y > as.numeric(U2[2]) - 0.015*4 & Y < as.numeric(U2[2]) + 0.015*4)

WUDG88_M2_dens <- M2_dens %>%
  filter(X > as.numeric(U1[1]) - 0.018*4 & X < as.numeric(U1[1]) + 0.018*4 &
           Y > as.numeric(U1[2]) - 0.018*4 & Y < as.numeric(U1[2]) + 0.018*4)

WUDG89_M2_dens <- M2_dens %>%
  filter(X > as.numeric(U2[1]) - 0.015*4 & X < as.numeric(U2[1]) + 0.015*4 &
           Y > as.numeric(U2[2]) - 0.015*4 & Y < as.numeric(U2[2]) + 0.015*4)

WUDG88_dens <- bind_rows(WUDG88_M1_dens, WUDG88_M2_dens) %>%
  mutate(Model = c(rep('Model 1', nrow(WUDG88_M1_dens)), rep('Model 2', nrow(WUDG88_M2_dens))))

tikz(file = "v11acs_WUD88_res.tex", standAlone=T, width = 6, height = 3)
ggplot() + geom_raster(data = WUDG88_dens, aes(X, Y, fill = z)) + 
  geom_sf(data = HPD_dens, aes(color = prob), fill = NA, size = 0.5) +
  coord_sf() +  theme_minimal() + 
  xlim(c(as.numeric(U1[1]) - 0.018*4,as.numeric(U1[1]) + 0.018*4)) + 
  ylim(c(as.numeric(U1[2]) - 0.018*4,as.numeric(U1[2]) + 0.018*4)) +
  facet_grid(~Model) + scale_fill_viridis_c(option = 'D', name = 'Density',
                                            end = max(WUDG88_dens$z)/max(dens$z)) + 
  scale_color_viridis_c(option = 'A', direction = -1, name = '$1-\\alpha$') +
  annotate('path', 
           x = as.numeric(U1[1]) + 0.018*cos(seq(0, 2*pi, length.out = 100)),
           y = as.numeric(U1[2]) + 0.018*sin(seq(0, 2*pi, length.out = 100)), 
           color = 'gold', size = 0.5, linetype = 'dotted') +
  xlab('X') + ylab('Y')
dev.off()
system('pdflatex v11acs_WUD88_res.tex')

WUDG89_dens <- bind_rows(WUDG89_M1_dens, WUDG89_M2_dens) %>%
  mutate(Model = c(rep('Model 1', nrow(WUDG89_M1_dens)), rep('Model 2', nrow(WUDG89_M2_dens))))

tikz(file = "v11acs_WUD89_res.tex", standAlone=T, width = 6, height = 3)
ggplot() + geom_raster(data = WUDG89_dens, aes(X, Y, fill = z)) + 
  geom_sf(data = HPD_dens, aes(color = prob), fill = NA, size = 0.5) +
  coord_sf() +  theme_minimal() + 
  xlim(c(as.numeric(U2[1]) - 0.015*4,as.numeric(U2[1]) + 0.015*4)) + 
  ylim(c(as.numeric(U2[2]) - 0.015*4,as.numeric(U2[2]) + 0.015*4)) +
  facet_grid(~Model) + scale_fill_viridis_c(option = 'D', name = '  Density',
                                            end = max(WUDG89_dens$z)/max(dens$z)) + 
  scale_color_viridis_c(option = 'A', direction = -1, name = '$1-\\alpha$') +
  annotate('path', 
           x = as.numeric(U2[1]) + 0.015*cos(seq(0, 2*pi, length.out = 100)),
           y = as.numeric(U2[2]) + 0.015*sin(seq(0, 2*pi, length.out = 100)), 
           color = 'gold', size = 0.5, linetype = 'dotted') +
  xlab('X') + ylab('Y')
dev.off()
system('pdflatex v11acs_WUD89_res.tex')




