library(tidyverse)
library(modeest)
library(coda)
library(stableGR)
library(RColorBrewer)
library(reshape2)
library(emdbook)
library(swfscMisc)
library(tikzDevice)
library(expint)
library(gplots)

options(warn = -1)
source('src/help_functions.R')
v_acs <- read_csv('data/v7acs_GC.csv')
col <- readRDS('data/v7_acs_col.RDS')
pts <- v_acs[,3:4]/4300

# chain for model 1
res <- readRDS('data/v7acs_res/v7acs-pcp-2022-08-03-seed-538665-100k.RDS')
# chain for model 2
res_col <- readRDS('data/v7acs_res/v7acs-pcp-col-2022-08-03-seed-538665-100k.RDS')

# compute membership probabilities
# extract parameters to compute P(Z != 0 | X) (first 10000 as bur-in)
sub_res <- res[[1]]
sub_fix <- res[[2]]

b0 <- exp(sub_fix[,'b0'])

prob <- matrix(0, 61, 100001)
for (j in 1:100001) {
  prob[,j] <- PZ_X(X = pts, UDG_par = sub_res[[j]], b0 = b0[j]) 
}
prob_ncol <- t(prob)

# extract parameters to compute P(Z |!= 0 | X,V) (first 10000 as bur-in)
sub_res_col <- res_col[[1]]
sub_fix_col <- res_col[[2]]

b0_col <- exp(sub_fix_col[,'b0'])
s0_col <- exp(sub_fix_col[,'sigma0'])

prob <- matrix(0, 61, 100001)
for (j in 1:100001) {
  prob[,j] <- PZ_XY(X = pts, Y = col, UDG_par = sub_res_col[[j]], b0 = b0_col[j], sigma0 = s0_col[j]) 
}
prob_col <- t(prob)

prob_ncol <- prob_ncol[-c(1:10000),]
prob_col <- prob_col[-c(1:10000),]

PSRF_ncol <- stable.GR(prob_ncol, multivariate = F)
PSRF_col <- stable.GR(prob_col, multivariate = F)

# posterior density plots
ncol_mod <- as.data.frame(exp(res[[2]][-c(1:10000),]))
prior_df_ncol <- data.frame(value = c(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 
                                      seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100)),
                            density = c(dlnorm(seq(min(ncol_mod$b0), max(ncol_mod$b0), length= 100), 3.9772319, 0.25),
                                        dlnorm(seq(min(ncol_mod$lc), max(ncol_mod$lc), length= 100), 0, 0.2)),
                            key = c(rep('b0', 100), rep('lc', 100)))

col_mod <- as.data.frame(exp(res_col[[2]][-c(1:10000),]))
prior_df_col <- data.frame(value = c(seq(min(col_mod$b0), max(col_mod$b0), length= 100), 
                                     seq(min(col_mod$lc), max(col_mod$lc), length= 100),
                                     seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100)),
                           density = c(dlnorm(seq(min(col_mod$b0), max(col_mod$b0), length= 100), 3.9772319, 0.25),
                                       dlnorm(seq(min(col_mod$lc), max(col_mod$lc), length= 100), 0, 0.2),
                                       dgamma(seq(min(col_mod$sigma0), max(col_mod$sigma0), length= 100), 0.2, 0.05)),
                           key = c(rep('b0', 100), rep('lc', 100), rep('sigma0', 100)))

ncol.labs <- c('$\\beta_0$', '$\\lambda_c$')
names(ncol.labs) <- c('b0', 'lc')

col.labs <- c('$\\beta_0$', '$\\lambda_c$', '$\\sigma_0$')
names(col.labs) <- c('b0', 'lc', 'sigma0')

tikz(file = "v7acs_post_ncol.tex", standAlone=T, width = 6, height = 3)
ggplot(gather(ncol_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_ncol, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol.labs))
dev.off()
system('pdflatex v7acs_post_ncol.tex')

tikz(file = "v7acs_post_col.tex", standAlone=T, width = 9, height = 3)
ggplot(gather(col_mod), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 50) + 
  geom_line(data = prior_df_col, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col.labs))
dev.off()
system('pdflatex v7acs_post_col.tex')

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

set.seed(781025)
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

U1 <- data.frame(x = 661, y = 2987)/4300
U2 <- data.frame(x = 1756, y = 752)/4300

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

tikz(file = "v7acs_U1_post_ncol.tex", standAlone=T, width = 9, height = 3)
ggplot(gather(df_U1_ncol), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_ncol_U1, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol_U1_labs))
dev.off()
system('pdflatex v7acs_U1_post_ncol.tex')

prior_df_ncol_U2 <- data.frame(value = c(seq(min(df_U2_ncol$N), max(df_U2_ncol$N), length= 100), 
                                         seq(min(df_U2_ncol$R), max(df_U2_ncol$R), length= 100),
                                         seq(min(df_U2_ncol$n), max(df_U2_ncol$n), length= 100)),
                               density = c(dlnorm(seq(min(df_U2_ncol$N), max(df_U2_ncol$N), length= 100), log(6), 1),
                                           dlnorm(seq(min(df_U2_ncol$R), max(df_U2_ncol$R), length= 100), -3.9, 0.5),
                                           dlnorm(seq(min(df_U2_ncol$n), max(df_U2_ncol$n), length= 100), 0, 0.75)),
                               key = c(rep('N', 100), rep('R', 100), rep('n', 100)))

ncol_U2_labs <- c('$N_u^2$', '$R_u^2$', '$n_u^2$')
names(ncol_U2_labs) <- c('N', 'R', 'n')

tikz(file = "v7acs_U2_post_ncol.tex", standAlone=T, width = 9, height = 3)
ggplot(gather(df_U2_ncol), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) +  
  geom_line(data = prior_df_ncol_U2, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = ncol_U2_labs))
dev.off()
system('pdflatex v7acs_U2_post_ncol.tex')

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

tikz(file = "v7acs_U1_post_col.tex", standAlone=T, width = 6, height = 6)
ggplot(gather(df_U1_col), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_col_U1, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col_U1_labs))
dev.off()
system('pdflatex v7acs_U1_post_col.tex')

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

tikz(file = "v7acs_U2_post_col.tex", standAlone=T, width = 6, height = 6)
ggplot(gather(df_U2_col), aes(value)) + 
  geom_histogram(aes(y = after_stat(density)), bins = 20) + 
  geom_line(data = prior_df_col_U2, aes(value, density), color = 'red') + 
  facet_wrap(~key, scales = 'free', labeller = labeller(key = col_U2_labs))
dev.off()
system('pdflatex v7acs_U2_post_col.tex')

# density and probability plots
ncol_dat <- res[[1]][-c(1:10000)]
col_dat <- res_col[[1]][-c(1:10000)]

C_post_r <- lapply(ncol_dat, function(x) x)[seq(1, 90001, by = 5)]
X_nc <- do.call(rbind, C_post_r)
X_nc <- as.data.frame(X_nc)

C_post_rc <- lapply(col_dat, function(x) x)[seq(1, 90001, by = 5)]
X_c <- do.call(rbind, C_post_rc)
X_c <- as.data.frame(X_c)

U <- rbind(U1, U2)

ci <- c(0.1, 0.2, 0.3, 0.4)
dm <- sapply(ci, function(x) detect_metric(U, c(0.05/3, 0.008), x, X_nc, X_c))
acc_bs <- vector('list', length = 4)
prec_bs <- vector('list', length = 4)
FP_bs <- vector('list', length = 4)
set.seed(690537)
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
    pm <- detect_metric(U, c(0.05/3, 0.008), ci[i], X_nc_j, X_c_j)
    acc_bs[[i]][j] <- pm[[1]]
    FP_bs[[i]][j] <- pm[[2]]
    prec_bs[[i]][j] <- pm[[3]]
  }
}

saveRDS(list(dm, acc_bs, FP_bs, prec_bs), 'data/v7acs_res/v7acs_pm.RDS')

dat_dm <- data.frame(prob = 1 - ci, acc = unlist(dm[1,]), FP = unlist(dm[2,]), prec = unlist(dm[3,]), 
                     TP_upper = unlist(lapply(acc_bs, function(x) quantile(x, 0.975))), 
                     TP_lower = unlist(lapply(acc_bs, function(x) quantile(x, 0.025))),
                     FP_upper = unlist(lapply(FP_bs, function(x) quantile(x, 0.975))), 
                     FP_lower = unlist(lapply(FP_bs, function(x) quantile(x, 0.025))),
                     prec_upper = unlist(lapply(prec_bs, function(x) quantile(x, 0.975))),
                     prec_lower = unlist(lapply(prec_bs, function(x) quantile(x, 0.025))))

ggplot(dat_dm, aes(prob, prec)) + geom_line(size = 0.2, color = "#A25050") + geom_point(size = 0.2, color = "#7C0000") +
  geom_errorbar(size = 0.1, aes(ymin = prec_lower, ymax = prec_upper), color = "#7C0000", width = 0.005) + 
  ylim(c(0.8,1.5)) + geom_abline(intercept = 1, slope = 0, linetype = 'dashed', size = 0.2) + theme_minimal()

ggplot(C_post_rc, aes(x, y)) +  geom_bin_2d(binwidth = 0.0025) + 
  scale_fill_viridis_c(option = 'D', name = 'Density') + xlim(c(0,1)) + ylim(c(0,1)) 

options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))

tikz(file = "v7acs_res.tex", standAlone=T, width = 7, height = 3.5)
ggplot(Xc, aes(x, y)) + stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE, n = 200) + 
  facet_wrap(~Model) +
  scale_fill_viridis_c(option = 'D', name = 'Density') + xlim(c(0,1)) + ylim(c(0,1)) + 
  geom_point(data = prob, aes(x,y, color = cut(value, c(0, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9,1))), size = 0.5) + facet_wrap(~Model) +
  scale_color_viridis_d(option = 'C', name = '$\\mathbb{P}(Z \\neq 0 \\mid \\mathrm{Data})$', direction = -1) + coord_fixed() + theme_minimal() +
  geom_point(data = U1, aes(x,y), size = 9, shape = 1, color = 'gold') + 
  geom_point(data = U2, aes(x,y), size = 6, shape = 1, color = 'gold') + 
  xlab('X') + ylab('Y')
dev.off()
system('pdflatex v7acs_res.tex')



