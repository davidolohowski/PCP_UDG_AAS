library(tidyverse)
library(mcmcse)
library(coda)
library(stableGR)
library(RColorBrewer)
library(reshape2)
library(emdbook)
library(swfscMisc)
library(tikzDevice)
library(gplots)

options(warn = -1)
source('src/help_functions.R')
v_acs <- read_csv('data/v13acs_GC.csv')
col <- readRDS('data/v13_acs_col.RDS')
pts <- v_acs[,3:4]/4300

# chain for model 1
res <- readRDS('data/v13acs_res/v13acs-pcp-2022-07-22-seed-409635-100k.RDS')
# chain for model 2
res_col <- readRDS('data/v13acs_res/v13acs-pcp-col-2022-07-22-seed-409635-100k.RDS')

# compute membership probabilities
# extract parameters to compute P(Z != 0 | X) (first 10000 as bur-in)
sub_res <- res[[1]]
sub_fix <- res[[2]]

b0 <- exp(sub_fix[,'b0'])

prob <- matrix(0, 32, 100001)
for (j in 1:100001) {
  prob[,j] <- PZ_X(X = pts, UDG_par = sub_res[[j]], b0 = b0[j]) 
}
prob_ncol <- t(prob)

# extract parameters to compute P(Z |!= 0 | X,V) (first 10000 as bur-in)
sub_res_col <- res_col[[1]]
sub_fix_col <- res_col[[2]]

b0_col <- exp(sub_fix_col[,'b0'])
s0_col <- exp(sub_fix_col[,'sigma0'])

prob <- matrix(0, 32, 100001)
for (j in 1:100001) {
  prob[,j] <- PZ_XY(X = pts, Y = col, UDG_par = sub_res_col[[j]], b0 = b0_col[j], sigma0 = s0_col[j]) 
}
prob_col <- t(prob)

prob_ncol <- prob_ncol[-c(1:10000),]
prob_col <- prob_col[-c(1:10000),]

PSRF_ncol <- stable.GR(prob_ncol, multivariate = F)
PSRF_col <- stable.GR(prob_col, multivariate = F)

ncol_mod <- as.data.frame(exp(res[[2]][-c(1:10000),]))
prior_df_ncol <- data.frame(value = c(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 
                                      seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100),
                                      seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100),
                                      seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100),
                                      seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100)),
                            density = c(dlnorm(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 4.487540, 0.25),
                                        dlnorm(seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100), log(1), 0.2),
                                        dlnorm(seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100), log(0.5), 0.5),
                                        dlnorm(seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100), log(20), 0.25),
                                        dlnorm(seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100), 4.6 - log(4300), 0.25)),
                            key = c(rep('b0', 100), rep('lc', 100), rep('n1', 100), rep('N1', 100),
                                    rep('R1', 100)))

col_mod <- as.data.frame(exp(res_col[[2]][-c(1:10000),]))
prior_df_col <- data.frame(value = c(seq(min(col_mod$b0), max(col_mod$b0), length= 100), 
                                     seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100),
                                     seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100),
                                     seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100),
                                     seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100),
                                     seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100)),
                           density = c(dlnorm(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 4.487540, 0.25),
                                       dlnorm(seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100), log(1), 0.2),
                                       dlnorm(seq(min(ncol_mod$n1), max(ncol_mod$n1), length= 100), log(0.5), 0.5),
                                       dlnorm(seq(min(ncol_mod$N1), max(ncol_mod$N1), length= 100), log(20), 0.25),
                                       dlnorm(seq(min(ncol_mod$R1), max(ncol_mod$R1), length= 100), 4.6 - log(4300), 0.25),
                                       dgamma(seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100), 0.2, 0.05)),
                           key = c(rep('b0', 100), rep('lc', 100), rep('n1', 100), rep('N1', 100),
                                   rep('R1', 100), rep('sigma0', 100)))

ncol.labs <- c('$\\beta_0$', '$\\lambda_c$', '$N_g^1$', '$n_g^1$', '$R_g^1$')
names(ncol.labs) <- c('b0', 'lc', 'N1', 'n1', 'R1')

col.labs <- c('$\\beta_0$', '$\\lambda_c$', '$N_g^1$', '$n_g^1$', '$R_g^1$', '$\\sigma_0$')
names(col.labs) <- c('b0', 'lc', 'N1', 'n1', 'R1', 'sigma0')

tikz(file = "v13acs_post_ncol.tex", standAlone=T, width = 9, height = 4.5)
ggplot(gather(ncol_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_ncol, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol.labs))
dev.off()
system('pdflatex v13acs_post_ncol.tex')

tikz(file = "v13acs_post_col.tex", standAlone=T, width = 9, height = 4.5)
ggplot(gather(col_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_col, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col.labs))
dev.off()
system('pdflatex v13acs_post_col.tex')


# density and probability plots
ncol_dat <- res[[1]][-c(1:10000)]
col_dat <- res_col[[1]][-c(1:10000)]

C_post_r <- lapply(ncol_dat, function(x) x)
C_post_r <- do.call(rbind, C_post_r)
C_post_r <- as.data.frame(C_post_r)

C_post_rc <- lapply(col_dat, function(x) x)
C_post_rc <- do.call(rbind, C_post_rc)
C_post_rc <- as.data.frame(C_post_rc)

N_nc <- 1*(unlist(lapply(ncol_dat, nrow)) == 0)
N_c <- 1*(unlist(lapply(col_dat, nrow)) == 0)

sum(N_nc)/90001
sum(N_c)/90001

p0_nc <- numeric(100)
p0_c <- numeric(100)

pu_nc <- numeric(100)
pu_c <- numeric(100)

hist_res_nc <- hist2d(C_post_r$x, C_post_r$y, nbins = c(7,7))
hist_res_c <- hist2d(C_post_rc$x, C_post_rc$y, nbins = c(7,7))

max(hist_res_nc$counts)/nrow(C_post_r)*(1 - sum(N_nc)/90001)
max(hist_res_c$counts)/nrow(C_post_rc)*(1 - sum(N_c)/90001)

set.seed(605233)
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
