library(tikzDevice)
library(tidyverse)

no_UDG <- data.frame(est = c(0.2661, 0.2664, 0.2325, 0.2683, 0.2611, 0.2745, 0.3317, 0.2941, 0.3404, 0.3487,
                             0, 0, 0.0278, 0.0503, 0.0001, 0, 0.0262, 0.0007, 0.0157, 0.0828,
                             0.0699, 0.0858, 0.1273, 0.1343, 0.0509, 0.0630, 0.0302, 0.0536, 0.0458, 0.0320,
                             0.4358, 0.4123, 0.4141, 0.4054, 0.4089, 0.6618, 0.1943, 0.2843, 0.5001, 0.4820),
                     lower = c(0.2564, 0.2566, 0.2225, 0.2563, 0.2501, 0.2636, 0.3246, 0.2862, 0.3295, 0.3402,
                               0, 0, 0.0211, 0.0435, 0, 0, 0.0219, 0.0002, 0.0099, 0.0703,
                               0.0609, 0.0741, 0.1159, 0.1217, 0.0428, 0.055, 0.0254, 0.0425, 0.0381, 0.0269,
                               0.4288, 0.3681, 0.3963, 0.3901, 0.2965, 0.5307, 0.1827, 0.2758, 0.4241, 0.4595),
                     upper = c(0.2751, 0.2763, 0.2413, 0.2829, 0.2706, 0.2822, 0.3384, 0.3057, 0.3520, 0.3564,
                               0, 0, 0.0352, 0.0599, 0.0004, 0, 0.0304, 0.0013, 0.0238, 0.0957,
                               0.0766, 0.0943, 0.1453, 0.1793, 0.0597, 0.0746, 0.0339, 0.0605, 0.0562, 0.0368,
                               0.4458, 0.4159, 0.4305, 0.4239, 0.4429, 0.6763, 0.2084, 0.2946, 0.5117, 0.5006),
                     pointings = rep(rep(c('V6-WFC3', 'V11-WFC3', 'V13-ACS', 'V13-WFC3', '14-WFC3', 'V7-ACS', 'V8-WFC3', 'V10-ACS', 'V11-ACS', 'V15-ACS'), each = 2), 2),
                     prob_type = rep(c('$\\hat{\\mathbb{P}}(N_c = 0 \\mid \\mathbf{D})$', '$\\max_i\\hat{\\mathbb{P}}(\\mathrm{UDGs} \\in \\omega_i \\mid \\mathbf{D})$'), each = 20),
                     UDG = rep(rep(c('Pointing has known UDGs? No', 'Pointing has known UDGs? Yes'), each = 10), 2),
                     Model = rep(rep(c('Model 1', 'Model 2'), 10), 2))

no_UDG$pointings <- factor(no_UDG$pointings, levels = rev(c('V6-WFC3', 'V11-WFC3', 'V13-ACS', 'V13-WFC3', '14-WFC3', 'V7-ACS', 'V8-WFC3', 'V10-ACS', 'V11-ACS', 'V15-ACS')))

options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ),"\\usepackage{amsfonts}"))

tikz(file = "UDG_existence.tex", standAlone=T, width = 9/1.2, height = 3/1.2)
ggplot(no_UDG, aes(est, pointings,color = prob_type, shape = Model)) + 
  geom_point(position = position_dodge(width = 0.5), size = 0.75) +
  geom_errorbar(aes(xmin = lower, xmax = upper), size = 0.5, width = 0.15, position = position_dodge(width = 0.5)) +
  facet_wrap(.~UDG, scales = 'free', drop = T) + theme_minimal() + 
  labs(x = 'Estimates', y = 'Pointing ID') + theme(legend.title=element_blank())
dev.off()
system('pdflatex UDG_existence.tex')
