library(tidyverse)
library(modeest)
library(coda)
library(stableGR)
library(RColorBrewer)
library(reshape2)
library(emdbook)
library(swfscMisc)
library(tikzDevice)
library(plotly)

options(warn = -1)
source('src/help_functions.R')
v_acs <- read_csv('data/v14acs_GC.csv')
col <- readRDS('data/v14_acs_col.RDS')
pts <- v_acs[,3:4]/4300

# chain for model 1
res <- readRDS('data/v14acs_res/v14acs-pcp-2022-11-30-seed-709568-100k.RDS')

# chain for model 2
res_col <- readRDS('data/v14acs_res/v14acs-pcp-col-2022-11-30-seed-709568-100k.RDS')

gal_fix <- data.frame(x = c(3368)/4300, y = c(2716)/4300,
                      theta = c(-10*pi/180), e = c(1.25))

# extract parameters to compute P(Z != 0 | X) (first 10000 as bur-in)
sub_res <- res[[1]]
sub_fix <- res[[2]]

b0 <- exp(sub_fix[,'b0'])
gal1_par <- exp(sub_fix[,3:5])

colnames(gal1_par) <- c('N', 'R', 'n')

gal1_par <- as.data.frame(gal1_par)

gal1_par$x <- gal_fix$x[1]
gal1_par$y <- gal_fix$y[1]
gal1_par$theta <- gal_fix$theta[1]
gal1_par$e <- gal_fix$e[1]

gal_par <- list()

for(j in 1:100001){
  gal_par[[j]] <- rbind(gal1_par[j,])
}

prob <- matrix(0, 182, 100001)
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

colnames(gal1_par_col) <- c('N', 'R', 'n')

gal1_par_col <- as.data.frame(gal1_par_col)

gal1_par_col$x <- gal_fix$x[1]
gal1_par_col$y <- gal_fix$y[1]
gal1_par_col$theta <- gal_fix$theta[1]
gal1_par_col$e <- gal_fix$e[1]

gal_par <- list()

for(j in 1:100001){
  gal_par[[j]] <- rbind(gal1_par_col[j,])
}

prob <- matrix(0, 182, 100001)
for (j in 1:100001) {
  prob[,j] <- PZ_XY(pts, col, sub_res_col[[j]], gal_par[[j]], b0_col[j], s0_col[j]) 
}

prob_col <- t(prob)

prob_ncol <- prob_ncol[-c(1:10001),]
prob_col <- prob_col[-c(1:10001),]

PSRF_ncol <- stable.GR(prob_ncol, multivariate = F)
PSRF_col <- stable.GR(prob_col, multivariate = F)

# posterior density plots
ncol_mod <- as.data.frame(exp(res[[2]][-c(1:10000),]))
prior_df_ncol <- data.frame(value = c(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 
                                      seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100),
                                      seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100),
                                      seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100),
                                      seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100)),
                            density = c(dlnorm(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 4.017537, 0.25),
                                        dlnorm(seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100), log(2), 0.2),
                                        dlnorm(seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100), log(0.5), 0.5),
                                        dlnorm(seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100), log(130), 0.25),
                                        dlnorm(seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100), 6.5 - log(4300), 0.25)),
                            key = c(rep('b0', 100), rep('lc', 100), rep('n1', 100), rep('N1', 100),
                                    rep('R1', 100)))

col_mod <- as.data.frame(exp(res_col[[2]][-c(1:10000),]))
prior_df_col <- data.frame(value = c(seq(min(col_mod$b0), max(col_mod$b0), length= 100), 
                                     seq(min(col_mod$lc), max(col_mod$lc), length= 100),
                                     seq(min(col_mod$n1), max(col_mod$n1), length= 100),
                                     seq(min(col_mod$N1), max(col_mod$N1), length= 100),
                                     seq(min(col_mod$R1), max(col_mod$R1), length= 100),
                                     seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100)),
                           density = c(dlnorm(seq(min(col_mod$b0), max(col_mod$b0), length= 100), 4.017537, 0.25),
                                       dlnorm(seq(min(col_mod$lc), max(col_mod$lc), length= 100), log(2), 0.2),
                                       dlnorm(seq(min(col_mod$n1), max(col_mod$n1), length= 100), log(0.5), 0.5),
                                       dlnorm(seq(min(col_mod$N1), max(col_mod$N1), length= 100), log(130), 0.25),
                                       dlnorm(seq(min(col_mod$R1), max(col_mod$R1), length= 100), 6.5 - log(4300), 0.25),
                                       dgamma(seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100), 0.2, 0.05)),
                           key = c(rep('b0', 100), rep('lc', 100), rep('n1', 100), rep('N1', 100),
                                   rep('R1', 100), rep('sigma0', 100)))

ncol.labs <- c('$\\beta_0$', '$\\lambda_c$', '$n^1_g$', '$N_g^1$', '$R_g^1$')
names(ncol.labs) <- c('b0', 'lc', 'n1', 'N1', 'R1')

col.labs <- c('$\\beta_0$', '$\\lambda_c$', '$n^1_g$', '$N_g^1$', '$R_g^1$', '$\\sigma_0$')
names(col.labs) <- c('b0', 'lc', 'n1', 'N1', 'R1', 'sigma0')

tikz(file = "v14acs_post_ncol.tex", standAlone=T, width = 9, height = 9)
ggplot(gather(ncol_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_ncol, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol.labs))
dev.off()
system('pdflatex v14acs_post_ncol.tex')

tikz(file = "v14acs_post_col.tex", standAlone=T, width = 9, height = 9)
ggplot(gather(col_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_col, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col.labs))
dev.off()
system('pdflatex v14acs_post_col.tex')

U <- data.frame(x = 1087, y = 918)/4300

df_U_ncol <- do.call(rbind, res[[1]][-c(1:10001)]) %>%
  filter((x - U$x)^2 + (y - U$y)^2 <= 0.01^2) %>%
  dplyr::select(N, R, n)

df_U_col <- do.call(rbind, res_col[[1]][-c(1:10001)]) %>%
  filter((x - U$x)^2 + (y - U$y)^2 <= 0.01^2)%>%
  dplyr::select(N, R, n, sigma)

prior_df_ncol_U <- data.frame(value = c(seq(min(df_U_ncol$N), max(df_U_ncol$N), length= 100), 
                                        seq(min(df_U_ncol$R), max(df_U_ncol$R), length= 100),
                                        seq(min(df_U_ncol$n), max(df_U_ncol$n), length= 100)),
                              density = c(dlnorm(seq(min(df_U_ncol$N), max(df_U_ncol$N), length= 100), log(3), 1),
                                          dlnorm(seq(min(df_U_ncol$R), max(df_U_ncol$R), length= 100), -3.9, 0.5),
                                          dlnorm(seq(min(df_U_ncol$n), max(df_U_ncol$n), length= 100), 0, 0.75)),
                              key = c(rep('N', 100), rep('R', 100), rep('n', 100)))

ncol_U_labs <- c('$N_u^1$', '$R_u^1$', '$n_u^1$')
names(ncol_U_labs) <- c('N', 'R', 'n')

tikz(file = "v14acs_U_post_ncol.tex", standAlone=T, width = 9, height = 3)
ggplot(gather(df_U_ncol), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_ncol_U, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol_U_labs))
dev.off()
system('pdflatex v14acs_U_post_ncol.tex')

prior_df_col_U <- data.frame(value = c(seq(min(df_U_col$N), max(df_U_col$N), length= 100), 
                                       seq(min(df_U_col$R), max(df_U_col$R), length= 100),
                                       seq(min(df_U_col$n), max(df_U_col$n), length= 100),
                                       seq(min(df_U_col$sigma), max(df_U_col$sigma), length= 100)),
                             density = c(dlnorm(seq(min(df_U_col$N), max(df_U_col$N), length= 100), log(3), 1),
                                         dlnorm(seq(min(df_U_col$R), max(df_U_col$R), length= 100), -3.9, 0.5),
                                         dlnorm(seq(min(df_U_col$n), max(df_U_col$n), length= 100), 0, 0.75),
                                         5/gamma(0.05)*gammainc(0.05 - 1, seq(min(df_U_col$sigma), max(df_U_col$sigma), length= 100))),
                             key = c(rep('N', 100), rep('R', 100), rep('n', 100), rep('sigma', 100)))

col_U_labs <- c('$N_u^1$', '$R_u^1$', '$n_u^1$', '$\\sigma_1$')
names(col_U_labs) <- c('N', 'R', 'n', 'sigma')

tikz(file = "v14acs_U_post_col.tex", standAlone=T, width = 6, height = 6)
ggplot(gather(df_U_col), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_col_U, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col_U_labs))
dev.off()
system('pdflatex v14acs_U_post_col.tex')


# density and probability plots
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

hist_res_nc <- hist2d(X_nc$x, X_nc$y, nbins = c(7,7))
hist_res_c <- hist2d(X_c$x, X_c$y, nbins = c(7,7))

max(hist_res_nc$counts)/nrow(X_nc)*(1 - sum(N_nc)/90001)
max(hist_res_c$counts)/nrow(X_c)*(1 - sum(N_c)/90001)

tikz(file = "v14acs_model1_res.tex", standAlone=T, width = 5, height = 4.5*5/6)
ggplot(X_nc, aes(x, y)) + stat_density_2d(aes(fill = ..density..), 
                                          geom = "raster", contour = FALSE, n = 200,
                                          h = c(bw.nrd(X_nc$x), bw.nrd(X_nc$y))) + 
  scale_fill_viridis_c(option = 'D', name = 'Density') + xlim(c(0,1)) + ylim(c(0,1)) + coord_fixed() +
  theme_minimal() + geom_point(data = pts, aes(x,y), size = 0.1)
dev.off()
system('pdflatex v14acs_model1_res.tex')

tikz(file = "v14acs_model2_res.tex", standAlone=T, width = 5, height = 4.5*5/6)
ggplot(X_c, aes(x, y)) + stat_density_2d(aes(fill = ..density..), 
                                         geom = "raster", contour = FALSE, n = 200,
                                         h = c(bw.nrd(X_c$x), bw.nrd(X_c$y))) + 
  scale_fill_viridis_c(option = 'D', name = 'Density') + xlim(c(0,1)) + ylim(c(0,1)) + coord_fixed()+
  theme_minimal() + geom_point(data = pts, aes(x,y), size = 0.1)
dev.off()
system('pdflatex v14acs_model2_res.tex')

set.seed(109857)
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
