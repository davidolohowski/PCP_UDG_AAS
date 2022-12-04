library(tidyverse)
library(tikzDevice)

f1 <- list()
f1_names <- list.files('data/PCP_sim_weak/')

f2 <- list()
f2_names <- list.files('data/PCP_sim_base/')

f3 <- list()
f3_names <- list.files('data/PCP_sim_strong/')

for (i in 1:100) {
  f1[[i]] <- readRDS(paste0('data/PCP_sim_weak/', f1_names[i]))
  f2[[i]] <- readRDS(paste0('data/PCP_sim_base/', f2_names[i]))
  f3[[i]] <- readRDS(paste0('data/PCP_sim_strong/', f3_names[i]))
}

prec1 <- (unlist(lapply(c(3, 8, 13, 18), function(y) {mean(unlist(lapply(f1, function(x) x[[y]])))})) - 1)*100
prec2 <- (unlist(lapply(c(3, 8, 13, 18), function(y) {mean(unlist(lapply(f2, function(x) x[[y]])))})) - 1)*100
prec3 <- (unlist(lapply(c(3, 8, 13, 18), function(y) {mean(unlist(lapply(f3, function(x) x[[y]])))})) - 1)*100

ps1 <- (unlist(lapply(c(3, 8, 13, 18), function(y) {sd(unlist(lapply(f1, function(x) x[[y]]*10)))})))
ps2 <- (unlist(lapply(c(3, 8, 13, 18), function(y) {sd(unlist(lapply(f2, function(x) x[[y]]*10)))})))
ps3 <- (unlist(lapply(c(3, 8, 13, 18), function(y) {sd(unlist(lapply(f3, function(x) x[[y]]*10)))})))

FP1 <- (unlist(lapply(c(2, 7, 12, 17), function(y) {mean(unlist(lapply(f1, function(x) x[[y]])))})))
FP2 <- (unlist(lapply(c(2, 7, 12, 17), function(y) {mean(unlist(lapply(f2, function(x) x[[y]])))})))
FP3 <- (unlist(lapply(c(2, 7, 12, 17), function(y) {mean(unlist(lapply(f3, function(x) x[[y]])))})))

fs1 <- (unlist(lapply(c(2, 7, 12, 17), function(y) {sd(unlist(lapply(f1, function(x) x[[y]])))})))/10
fs2 <- (unlist(lapply(c(2, 7, 12, 17), function(y) {sd(unlist(lapply(f2, function(x) x[[y]])))})))/10
fs3 <- (unlist(lapply(c(2, 7, 12, 17), function(y) {sd(unlist(lapply(f3, function(x) x[[y]])))})))/10

TP1 <- (unlist(lapply(c(1, 6, 11, 16), function(y) {mean(unlist(lapply(f1, function(x) x[[y]])))})))
TP2 <- (unlist(lapply(c(1, 6, 11, 16), function(y) {mean(unlist(lapply(f2, function(x) x[[y]])))})))
TP3 <- (unlist(lapply(c(1, 6, 11, 16), function(y) {mean(unlist(lapply(f3, function(x) x[[y]])))})))

ts1 <- (unlist(lapply(c(1, 6, 11, 16), function(y) {sd(unlist(lapply(f1, function(x) x[[y]])))})))/10
ts2 <- (unlist(lapply(c(1, 6, 11, 16), function(y) {sd(unlist(lapply(f2, function(x) x[[y]])))})))/10
ts3 <- (unlist(lapply(c(1, 6, 11, 16), function(y) {sd(unlist(lapply(f3, function(x) x[[y]])))})))/10


df <- data.frame(mean = c(TP1, TP2, TP3, prec1, prec2, prec3, FP1, FP2, FP3),
                 lower = c(TP1 - 1.96*ts1, TP2 - 1.96*ts2, TP3 - 1.96*ts3, 
                           prec1 - 1.96*ps1, prec2 - 1.96*ps2, prec3 - 1.96*ps3,
                           FP1 - 1.96*ts1, FP2 - 1.96*fs2, FP3 + 1.96*fs3),
                 upper = c(TP1 + 1.96*ts1, TP2 + 1.96*ts2, TP3 + 1.96*ts3, 
                           prec1 + 1.96*ps1, prec2 + 1.96*ps2, prec3 + 1.96*ps3,
                           FP1 + 1.96*ts1, FP2 + 1.96*fs2, FP3 + 1.96*fs3),
                 ss = rep(rep(c('Weak', 'Base', 'Strong'), each = 4),3),
                 alpha = rep(rep(as.character(c(0.1, 0.2, 0.3, 0.4)), 3), 3),
                 measure = rep(rep(c('TP Difference', 'Relative \\% Increase in Prec', 'FP Difference'), each = 4), each = 3))

df$ss <- factor(df$ss, levels = c('Weak', 'Base', 'Strong'))

options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))

tikz(file = "metric_compare_sim.tex", standAlone=T, width = 6*1.2, height = 2*1.2)
ggplot(df, aes(ss, mean)) + geom_line(aes(color = alpha, group = alpha), position = position_dodge(width = 0.25)) +
  geom_point(aes(color = alpha, group = alpha), position = position_dodge(width = 0.25), size = 0.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper, color = alpha, width = 0.05), position = position_dodge(width = 0.25)) + 
  theme_minimal() + 
  scale_color_viridis_d(name = '$\\alpha$') + facet_wrap(.~measure, scales = 'free_y') + labs(x = 'UDG signal strength', y = 'Model 1 vs. Model 2')
dev.off()
system('pdflatex metric_compare_sim.tex')
