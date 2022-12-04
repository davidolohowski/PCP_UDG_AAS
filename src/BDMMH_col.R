source('src/help_functions.R')

# Birth-death-move Metropolis Hasting algorithm for fitting Model 2 in the paper

BDM_MH_supperpose_color <- function(X, #point pattern
                                    Y, #colour variation
                                    S, #study region
                                    b0_start, #starting value for beta_0
                                    gal_fix = NULL, #fixed parameter values for normal galaxies
                                    gal_start = NULL, #starting values (log-scale) for unknown parameters associated with normal galaxies 
                                    UDG_start, #starting values (natural scale) for parameters of UDGs
                                    lc_start, #starting value (log-scale) for lambda_c
                                    s0_start, #starting value (log-scale) for sigma_0
                                    npix = 50000, # number of pixels that discretise the study region for integration
                                    gal_prior = NULL, # prior parameters for parameters of normal galaxies
                                    UDG_prior, # prior parameters for parameters of UDGs
                                    b0_prior, # prior parameters for beta_0
                                    lc_prior, # prior parameters for lambda_c
                                    s0_prior, # prior parameters for sigma_0
                                    M, # number of iterations
                                    norm_gal = FALSE){
  
  # construct integration grid for the integral in Poisson process likelihood
  grid <- sp::makegrid(S, n = npix)
  grid <- sp::SpatialPoints(grid, proj4string = CRS(proj4string(S)))
  grid <- raster::crop(grid, S)
  gridded(grid) <- T
  
  dx1 <- unname(grid@grid@cellsize[1])
  dx2 <- unname(grid@grid@cellsize[2])
  
  grid <- as.data.frame(grid)
  names(grid) <- c('x', 'y')
  grid <- as.matrix(grid)
  
  grid <- list(grid, c(dx1, dx2))
  
  # Model without normal galaxies
  if(!norm_gal){
    # initialize Theta
    fix_par <- matrix(c(b0_start, lc_start, s0_start), ncol = 3)
    colnames(fix_par) <- c('b0', 'lc', 'sigma0')
    
    # initialize with nc_start number of UDGs
    A <- sf::st_area(sf::st_as_sf(S))
    nc_start <- rpois(1, lc_start*A) 
    
    #initialize Phi
    phi <- list()
    if(nc_start == 0){
      phi[[1]] <- data.frame(x = numeric(0), y = numeric(0), 
                             N = numeric(0), R = numeric(0),
                             e = numeric(0), n = numeric(0), 
                             theta = numeric(0), sigma = numeric(0))
    }
    else{
      start_point <- as.data.frame(sp::spsample(S, nc_start, type = 'random'))
      
      phi[[1]] <- data.frame(x = start_point[,1], y = start_point[,2], 
                             N = rep(UDG_start$N, nc_start), R = rep(UDG_start$R, nc_start),
                             e = rep(UDG_start$e, nc_start), n = rep(UDG_start$n, nc_start), 
                             theta = rep(UDG_start$theta, nc_start), sigma = rep(UDG_start$sigma, nc_start))
    }
    
    # MCMC
    for (i in 1:M) {
      # Birth a new cluster if Phi is empty
      if(nrow(phi[[i]]) == 0){
        new_point <- as.data.frame(sp::spsample(S, 1, type = 'random'))
        new_phi <- data.frame(x = new_point[,1], y = new_point[,2], 
                              N = rlnorm(1, UDG_prior$N_mean, UDG_prior$N_sd), R = rlnorm(1, UDG_prior$R_mean, UDG_prior$R_sd),
                              e = rlnorm(1, UDG_prior$e_mean, UDG_prior$e_sd), n = rlnorm(1, UDG_prior$n_mean, UDG_prior$n_sd), 
                              theta = runif(1, 0, pi), sigma = runif(1, 0, exp(fix_par[i, 'sigma0'])))
        old_phi <- phi[[i]]
        # adpative MCMC for Theta
        if(i <= 1000){
          new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
          old_b <- fix_par[i, 'b0']
          new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
          old_lc <- fix_par[i, 'lc']
          new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
          old_s0 <- fix_par[i, 'sigma0']
        }
        else if(i <= 10000){
          new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par))
          new_b <- new_fix[1]
          old_b <- fix_par[i, 'b0']
          new_lc <- new_fix[2]
          old_lc <- fix_par[i, 'lc']
          new_s0 <- new_fix[3]
          old_s0 <- fix_par[i, 'sigma0']
        }
        else{
          new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par[1:10000,]))
          new_b <- new_fix[1]
          old_b <- fix_par[i, 'b0']
          new_lc <- new_fix[2]
          old_lc <- fix_par[i, 'lc']
          new_s0 <- new_fix[3]
          old_s0 <- fix_par[i, 'sigma0']
        }
        # compute MH ratio
        r_b <- log_post_superpose_color(X, Y, S, par_UDG = new_phi, sigma0 = new_s0, 
                                        b0 = new_b, lc = new_lc, grid = grid,
                                        UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                        lc_prior = lc_prior, s0_prior = s0_prior) + log(A) -
               log_post_superpose_color(X, Y, S, par_UDG = old_phi, sigma0 = old_s0, 
                                        b0 = old_b, lc = old_lc, grid = grid,
                                        UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                        lc_prior = lc_prior, s0_prior = s0_prior) - 
               dlnorm(new_phi$N, UDG_prior$N_mean, UDG_prior$N_sd, TRUE) - 
               dlnorm(new_phi$R, UDG_prior$R_mean, UDG_prior$R_sd, TRUE) -
               dlnorm(new_phi$n, UDG_prior$n_mean, UDG_prior$n_sd, TRUE) - 
               dlnorm(new_phi$e, UDG_prior$e_mean, UDG_prior$e_sd, TRUE) + 
               log(pi) + fix_par[i, 'sigma0']
        # accept or reject
        if(log(runif(1)) < r_b){
          phi[[i+1]] <- new_phi
          new_fix <- c(new_b, new_lc, new_s0)
          fix_par <- rbind(fix_par, new_fix)
        }
        else{
          phi[[i+1]] <- old_phi
          old_fix <- c(old_b, old_lc, old_s0)
          fix_par <- rbind(fix_par, old_fix)
        }
      }
      else{
        U_m <- runif(1)
        
        # do not move 
        if(U_m > 0.5){
          # birth step
          if(runif(1) < 0.5){
            # birth new cluster
            new_point <- as.data.frame(sp::spsample(S, 1, type = 'random'))
            new_clust <- data.frame(x = new_point[,1], y = new_point[,2], 
                                    N = rlnorm(1, UDG_prior$N_mean, UDG_prior$N_sd), R = rlnorm(1, UDG_prior$R_mean, UDG_prior$R_sd),
                                    e = rlnorm(1, UDG_prior$e_mean, UDG_prior$e_sd), n = rlnorm(1, UDG_prior$n_mean, UDG_prior$n_sd), 
                                    theta = runif(1, 0, pi), sigma = runif(1, 0, exp(fix_par[i, 'sigma0'])))
            old_phi <- phi[[i]]
            new_phi <- bind_rows(old_phi, new_clust)
            # adaptive MCMC for Theta
            if(i <= 1000){
              new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
              old_b <- fix_par[i, 'b0']
              new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
              old_lc <- fix_par[i, 'lc']
              new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
              old_s0 <- fix_par[i, 'sigma0']
            }
            else if(i <= 10000){
              new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par))
              new_b <- new_fix[1]
              old_b <- fix_par[i, 'b0']
              new_lc <- new_fix[2]
              old_lc <- fix_par[i, 'lc']
              new_s0 <- new_fix[3]
              old_s0 <- fix_par[i, 'sigma0']
            }
            else{
              new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par[1:10000,]))
              new_b <- new_fix[1]
              old_b <- fix_par[i, 'b0']
              new_lc <- new_fix[2]
              old_lc <- fix_par[i, 'lc']
              new_s0 <- new_fix[3]
              old_s0 <- fix_par[i, 'sigma0']
            }
            # MH ratio
            r_b <- log_post_superpose_color(X, Y, S, par_UDG = new_phi, sigma0 = new_s0, 
                                            b0 = new_b, lc = new_lc, grid = grid,
                                            UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                            lc_prior = lc_prior, s0_prior = s0_prior) + log(A) -
                   log_post_superpose_color(X, Y, S, par_UDG = old_phi, sigma0 = old_s0, 
                                            b0 = old_b, lc = old_lc, grid = grid,
                                            UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                            lc_prior = lc_prior, s0_prior = s0_prior) - log(nrow(old_phi) + 1)- 
                   dlnorm(new_clust$N, UDG_prior$N_mean, UDG_prior$N_sd, TRUE) - 
                   dlnorm(new_clust$R, UDG_prior$R_mean, UDG_prior$R_sd, TRUE) -
                   dlnorm(new_clust$n, UDG_prior$n_mean, UDG_prior$n_sd, TRUE) - 
                   dlnorm(new_clust$e, UDG_prior$e_mean, UDG_prior$e_sd, TRUE) + 
                   log(pi) + fix_par[i, 'sigma0']
            #accept or reject
            if(log(runif(1)) < r_b){
              phi[[i+1]] <- new_phi
              new_fix <- c(new_b, new_lc, new_s0)
              fix_par <- rbind(fix_par, new_fix)
            }
            else{
              phi[[i+1]] <- old_phi
              old_fix <- c(old_b, old_lc, old_s0)
              fix_par <- rbind(fix_par, old_fix)
            }
          }
          # death step
          else{
            old_phi <- phi[[i]]
            #choose a cluster to remove
            idx <- sample(1:nrow(old_phi), size = 1)
            new_phi <- old_phi[-idx,]
            # adaptive MCMC for Theta
            if(i <= 1000){
              new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
              old_b <- fix_par[i, 'b0']
              new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
              old_lc <- fix_par[i, 'lc']
              new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
              old_s0 <- fix_par[i, 'sigma0']
            }
            else if(i <= 10000){
              new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par))
              new_b <- new_fix[1]
              old_b <- fix_par[i, 'b0']
              new_lc <- new_fix[2]
              old_lc <- fix_par[i, 'lc']
              new_s0 <- new_fix[3]
              old_s0 <- fix_par[i, 'sigma0']
            }
            else{
              new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par[1:10000,]))
              new_b <- new_fix[1]
              old_b <- fix_par[i, 'b0']
              new_lc <- new_fix[2]
              old_lc <- fix_par[i, 'lc']
              new_s0 <- new_fix[3]
              old_s0 <- fix_par[i, 'sigma0']
            }
            # MH ratio
            r_d <- -1*(log_post_superpose_color(X, Y, S, par_UDG = old_phi, sigma0 = old_s0, 
                                                b0 = old_b, lc = old_lc, grid = grid,
                                                UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                                lc_prior = lc_prior, s0_prior = s0_prior) + log(A) - 
                       log_post_superpose_color(X, Y, S, par_UDG = new_phi, sigma0 = new_s0, 
                                                b0 = new_b, lc = new_lc, grid = grid,
                                                UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                                lc_prior = lc_prior, s0_prior = s0_prior) - log(nrow(old_phi) + 1)- 
                       dlnorm(old_phi$N[idx], UDG_prior$N_mean, UDG_prior$N_sd, TRUE) - 
                       dlnorm(old_phi$R[idx], UDG_prior$R_mean, UDG_prior$R_sd, TRUE) -
                       dlnorm(old_phi$n[idx], UDG_prior$n_mean, UDG_prior$n_sd, TRUE) - 
                       dlnorm(old_phi$e[idx], UDG_prior$e_mean, UDG_prior$e_sd, TRUE) +
                       log(pi) + fix_par[i, 'sigma0'])
            #accept or reject
            if(log(runif(1)) < r_d){
              phi[[i+1]] <- new_phi
              new_fix <- c(new_b, new_lc, new_s0)
              fix_par <- rbind(fix_par, new_fix)
            }
            else{
              phi[[i+1]] <- old_phi
              old_fix <- c(old_b, old_lc, old_s0)
              fix_par <- rbind(fix_par, old_fix)
            }
          }
        }
        # move step
        else{
          old_phi <- phi[[i]]
          # choose a cluster to move
          idx <- sample(1:nrow(old_phi), size = 1)
          new_clust <- data.frame(x = rnorm(1, old_phi$x[idx], 0.01), y = rnorm(1, old_phi$y[idx], 0.01),
                                  N = rlnorm(1, log(old_phi$N[idx]), 0.05), R = rlnorm(1, log(old_phi$R[idx]), 0.05),
                                  e = rlnorm(1, log(old_phi$e[idx]), 0.05), n = rlnorm(1, log(old_phi$n[idx]), 0.05), 
                                  theta = rnorm(1, old_phi$theta[idx], 0.05), sigma = rlnorm(1, log(old_phi$sigma[idx]), 0.05))
          new_phi <- bind_rows(old_phi[-idx,], new_clust)
          # adaptive MCMC for Theta
          if(i <= 1000){
            new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
            old_b <- fix_par[i, 'b0']
            new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
            old_lc <- fix_par[i, 'lc']
            new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
            old_s0 <- fix_par[i, 'sigma0']
          }
          else if(i <= 10000){
            new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par))
            new_b <- new_fix[1]
            old_b <- fix_par[i, 'b0']
            new_lc <- new_fix[2]
            old_lc <- fix_par[i, 'lc']
            new_s0 <- new_fix[3]
            old_s0 <- fix_par[i, 'sigma0']
          }
          else{
            new_fix <- MASS::mvrnorm(1, fix_par[i,], cov(fix_par[1:10000,]))
            new_b <- new_fix[1]
            old_b <- fix_par[i, 'b0']
            new_lc <- new_fix[2]
            old_lc <- fix_par[i, 'lc']
            new_s0 <- new_fix[3]
            old_s0 <- fix_par[i, 'sigma0']
          }
          # MH ratio
          r_m <- log_post_superpose_color(X, Y, S, par_UDG = new_phi, sigma0 = new_s0, 
                                          b0 = new_b, lc = new_lc, grid = grid,
                                          UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                          lc_prior = lc_prior, s0_prior = s0_prior) -
                 log_post_superpose_color(X, Y, S, par_UDG = old_phi, sigma0 = old_s0, 
                                          b0 = old_b, lc = old_lc, grid = grid,
                                          UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                          lc_prior = lc_prior, s0_prior = s0_prior) +
                 dlnorm(old_phi$N[idx], log(new_clust$N), 0.05, TRUE) - dlnorm(new_clust$N, log(old_phi$N[idx]), 0.05, TRUE) +
                 dlnorm(old_phi$R[idx], log(new_clust$R), 0.05, TRUE) - dlnorm(new_clust$R, log(old_phi$R[idx]), 0.05, TRUE) +
                 dlnorm(old_phi$n[idx], log(new_clust$n), 0.05, TRUE) - dlnorm(new_clust$n, log(old_phi$n[idx]), 0.05, TRUE) +
                 dlnorm(old_phi$e[idx], log(new_clust$e), 0.05, TRUE) - dlnorm(new_clust$e, log(old_phi$e[idx]), 0.05, TRUE) +
                 dlnorm(old_phi$sigma[idx], log(new_clust$sigma), 0.05, TRUE) - dlnorm(new_clust$sigma, log(old_phi$sigma[idx]), 0.05, TRUE)
          # accept or reject
          if(log(runif(1)) < r_m){
            phi[[i+1]] <- new_phi
            new_fix <- c(new_b, new_lc, new_s0)
            fix_par <- rbind(fix_par, new_fix)
          }
          else{
            phi[[i+1]] <- old_phi
            old_fix <- c(old_b, old_lc, old_s0)
            fix_par <- rbind(fix_par, old_fix)
          }
        }
      }
    }
    return(list(phi, fix_par))  
  }
  # model with normal galaxies
  else{
    if(any(c(is.null(gal_fix), is.null(gal_start), is.null(gal_prior)))){
      stop('Need to supply all parameters and priors associated with normal galaxies.')
    }
    else{
      # initialize Theta
      fix_par <- matrix(c(b0_start, lc_start, s0_start), ncol = 3)
      for (i in 1:nrow(gal_start)) {
        gal_i <- matrix(gal_start[i,], ncol = 3)
        fix_par <- cbind(fix_par, gal_i)
      }
      
      name <- c('b0', 'lc', 'sigma0')
      for (i in 1:nrow(gal_start)) {
        name <- c(name, c(paste0('N',i), paste0('R',i), paste0('n',i)))
      }
      
      colnames(fix_par) <- name
      
      # generate nc_start number of clusters
      A <- sf::st_area(sf::st_as_sf(S))
      nc_start <- rpois(1, exp(lc_start)*A) 
      
      # initialize Phi
      phi <- list()
      if(nc_start == 0){
        phi[[1]] <- data.frame(x = numeric(0), y = numeric(0), 
                               N = numeric(0), R = numeric(0),
                               e = numeric(0), n = numeric(0), 
                               theta = numeric(0), sigma = numeric(0))
      }
      else{
        start_point <- as.data.frame(sp::spsample(S, nc_start, type = 'random'))
        
        phi[[1]] <- data.frame(x = start_point[,1], y = start_point[,2], 
                               N = rep(UDG_start$N, nc_start), R = rep(UDG_start$R, nc_start),
                               e = rep(UDG_start$e, nc_start), n = rep(UDG_start$n, nc_start), 
                               theta = rep(UDG_start$theta, nc_start), sigma = rep(UDG_start$sigma, nc_start))
      }
      # MCMC
      for (i in 1:M) {
        # birth a new cluster if Phi is empty
        if(nrow(phi[[i]]) == 0){
          new_point <- as.data.frame(sp::spsample(S, 1, type = 'random'))
          new_phi <- data.frame(x = new_point[,1], y = new_point[,2], 
                                N = rlnorm(1, UDG_prior$N_mean, UDG_prior$N_sd), R = rlnorm(1, UDG_prior$R_mean, UDG_prior$R_sd),
                                e = rlnorm(1, UDG_prior$e_mean, UDG_prior$e_sd), n = rlnorm(1, UDG_prior$n_mean, UDG_prior$n_sd), 
                                theta = runif(1, 0, pi), sigma = runif(1, 0, exp(fix_par[i, 'sigma0'])))
          old_phi <- phi[[i]]
          # Adaptive MCMC for Theta
          if(i <= 1000){
            new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
            old_b <- fix_par[i, 'b0']
            new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
            old_lc <- fix_par[i, 'lc']
            new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
            old_s0 <- fix_par[i, 'sigma0']
            new_g <- rnorm(3*nrow(gal_start), fix_par[i, 4:(3+3*nrow(gal_start))], 0.025)
            old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
            new_fix <- c(new_b, new_lc, new_s0, new_g)
          }
          else if(i <= 10000){
            new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par))
            new_b <- new_fix[1]
            old_b <- fix_par[i, 'b0']
            new_lc <- new_fix[2]
            old_lc <- fix_par[i, 'lc']
            new_s0 <- new_fix[3]
            old_s0 <- fix_par[i, 'sigma0']
            new_g <- new_fix[4:(3+3*nrow(gal_start))]
            old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
          }
          else{
            new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par[1:10000,]))
            new_b <- new_fix[1]
            old_b <- fix_par[i, 'b0']
            new_lc <- new_fix[2]
            old_lc <- fix_par[i, 'lc']
            new_s0 <- new_fix[3]
            old_s0 <- fix_par[i, 'sigma0']
            new_g <- new_fix[4:(3+3*nrow(gal_start))]
            old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
          }
          
          new_gal <- matrix(new_g, nrow = nrow(gal_start), byrow = T)
          colnames(new_gal) <- c('N', 'R', 'n')
          old_gal <- matrix(old_g, nrow = nrow(gal_start), byrow = T)
          colnames(old_gal) <- c('N', 'R', 'n')
          
          # MH ratio
          r_b <- log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = new_gal, 
                                          par_UDG = new_phi, sigma0 = new_s0, b0 = new_b, 
                                          lc = new_lc, grid = grid, gal_prior = gal_prior, 
                                          UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                          lc_prior = lc_prior, s0_prior = s0_prior) + log(A) -
                 log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = old_gal, 
                                          par_UDG = old_phi, sigma0 = old_s0, b0 = old_b, 
                                          lc = old_lc, grid = grid, gal_prior = gal_prior, 
                                          UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                          lc_prior = lc_prior, s0_prior = s0_prior) - 
                 dlnorm(new_phi$N, UDG_prior$N_mean, UDG_prior$N_sd, TRUE) - 
                 dlnorm(new_phi$R, UDG_prior$R_mean, UDG_prior$R_sd, TRUE) -
                 dlnorm(new_phi$n, UDG_prior$n_mean, UDG_prior$n_sd, TRUE) - 
                 dlnorm(new_phi$e, UDG_prior$e_mean, UDG_prior$e_sd, TRUE) + 
                 log(pi) + fix_par[i, 'sigma0']
          # accept or reject
          if(log(runif(1)) < r_b){
            phi[[i+1]] <- new_phi
            fix_par <- rbind(fix_par, new_fix)
          }
          else{
            phi[[i+1]] <- old_phi
            fix_par <- rbind(fix_par, fix_par[i,])
          }
        }
        # if Phi is not empty:
        else{
          U_m <- runif(1)
          
          # birth a new cluster
          if(U_m > 0.5){
            if(runif(1) < 0.5){
              new_point <- as.data.frame(sp::spsample(S, 1, type = 'random'))
              new_clust <- data.frame(x = new_point[,1], y = new_point[,2], 
                                      N = rlnorm(1, UDG_prior$N_mean, UDG_prior$N_sd), R = rlnorm(1, UDG_prior$R_mean, UDG_prior$R_sd),
                                      e = rlnorm(1, UDG_prior$e_mean, UDG_prior$e_sd), n = rlnorm(1, UDG_prior$n_mean, UDG_prior$n_sd), 
                                      theta = runif(1, 0, pi), sigma = runif(1, 0, exp(fix_par[i, 'sigma0'])))
              old_phi <- phi[[i]]
              new_phi <- bind_rows(old_phi, new_clust)
              # adaptive MCMC for Theta
              if(i <= 1000){
                new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
                old_b <- fix_par[i, 'b0']
                new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
                old_lc <- fix_par[i, 'lc']
                new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
                old_s0 <- fix_par[i, 'sigma0']
                new_g <- rnorm(3*nrow(gal_start), fix_par[i, 4:(3+3*nrow(gal_start))], 0.025)
                old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
                new_fix <- c(new_b, new_lc, new_s0, new_g)
              }
              else if(i <= 10000){
                new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par))
                new_b <- new_fix[1]
                old_b <- fix_par[i, 'b0']
                new_lc <- new_fix[2]
                old_lc <- fix_par[i, 'lc']
                new_s0 <- new_fix[3]
                old_s0 <- fix_par[i, 'sigma0']
                new_g <- new_fix[4:(3+3*nrow(gal_start))]
                old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
              }
              else{
                new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par[1:10000,]))
                new_b <- new_fix[1]
                old_b <- fix_par[i, 'b0']
                new_lc <- new_fix[2]
                old_lc <- fix_par[i, 'lc']
                new_s0 <- new_fix[3]
                old_s0 <- fix_par[i, 'sigma0']
                new_g <- new_fix[4:(3+3*nrow(gal_start))]
                old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
              }
              
              new_gal <- matrix(new_g, nrow = nrow(gal_start), byrow = T)
              colnames(new_gal) <- c('N', 'R', 'n')
              old_gal <- matrix(old_g, nrow = nrow(gal_start), byrow = T)
              colnames(old_gal) <- c('N', 'R', 'n')
              
              # MH ratio
              r_b <- log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = new_gal, 
                                              par_UDG = new_phi, sigma0 = new_s0, b0 = new_b, 
                                              lc = new_lc, grid = grid, gal_prior = gal_prior, 
                                              UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                              lc_prior = lc_prior, s0_prior = s0_prior) + log(A) -
                     log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = old_gal, 
                                              par_UDG = old_phi, sigma0 = old_s0, b0 = old_b, 
                                              lc = old_lc, grid = grid, gal_prior = gal_prior, 
                                              UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                              lc_prior = lc_prior, s0_prior = s0_prior) - log(nrow(old_phi) + 1)- 
                     dlnorm(new_clust$N, UDG_prior$N_mean, UDG_prior$N_sd, TRUE) - 
                     dlnorm(new_clust$R, UDG_prior$R_mean, UDG_prior$R_sd, TRUE) -
                     dlnorm(new_clust$n, UDG_prior$n_mean, UDG_prior$n_sd, TRUE) - 
                     dlnorm(new_clust$e, UDG_prior$e_mean, UDG_prior$e_sd, TRUE) +
                     log(pi) + fix_par[i, 'sigma0']
              # accept or reject
              if(log(runif(1)) < r_b){
                phi[[i+1]] <- new_phi
                fix_par <- rbind(fix_par, new_fix)
              }
              else{
                phi[[i+1]] <- old_phi
                fix_par <- rbind(fix_par, fix_par[i,])
              }
            }
            # death step
            else{
              old_phi <- phi[[i]]
              # choose a cluster to remove
              idx <- sample(1:nrow(old_phi), size = 1)
              new_phi <- old_phi[-idx,]
              # adaptive MCMC for Theta
              if(i <= 1000){
                new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
                old_b <- fix_par[i, 'b0']
                new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
                old_lc <- fix_par[i, 'lc']
                new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
                old_s0 <- fix_par[i, 'sigma0']
                new_g <- rnorm(3*nrow(gal_start), fix_par[i, 4:(3+3*nrow(gal_start))], 0.025)
                old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
                new_fix <- c(new_b, new_lc, new_s0, new_g)
              }
              else if(i <= 10000){
                new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par))
                new_b <- new_fix[1]
                old_b <- fix_par[i, 'b0']
                new_lc <- new_fix[2]
                old_lc <- fix_par[i, 'lc']
                new_s0 <- new_fix[3]
                old_s0 <- fix_par[i, 'sigma0']
                new_g <- new_fix[4:(3+3*nrow(gal_start))]
                old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
              }
              else{
                new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par[1:10000,]))
                new_b <- new_fix[1]
                old_b <- fix_par[i, 'b0']
                new_lc <- new_fix[2]
                old_lc <- fix_par[i, 'lc']
                new_s0 <- new_fix[3]
                old_s0 <- fix_par[i, 'sigma0']
                new_g <- new_fix[4:(3+3*nrow(gal_start))]
                old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
              }
              
              new_gal <- matrix(new_g, nrow = nrow(gal_start), byrow = T)
              colnames(new_gal) <- c('N', 'R', 'n')
              old_gal <- matrix(old_g, nrow = nrow(gal_start), byrow = T)
              colnames(old_gal) <- c('N', 'R', 'n')
              
              # MH ratio
              r_d <- -1*(log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = old_gal, 
                                                  par_UDG = old_phi, sigma0 = old_s0, b0 = old_b, 
                                                  lc = old_lc, grid = grid, gal_prior = gal_prior, 
                                                  UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                                  lc_prior = lc_prior, s0_prior = s0_prior) + log(A) - 
                         log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = new_gal, 
                                                  par_UDG = new_phi, sigma0 = new_s0, b0 = new_b, 
                                                  lc = new_lc, grid = grid, gal_prior = gal_prior, 
                                                  UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                                  lc_prior = lc_prior, s0_prior = s0_prior) - log(nrow(old_phi) + 1)- 
                         dlnorm(old_phi$N[idx], UDG_prior$N_mean, UDG_prior$N_sd, TRUE) - 
                         dlnorm(old_phi$R[idx], UDG_prior$R_mean, UDG_prior$R_sd, TRUE) -
                         dlnorm(old_phi$n[idx], UDG_prior$n_mean, UDG_prior$n_sd, TRUE) - 
                         dlnorm(old_phi$e[idx], UDG_prior$e_mean, UDG_prior$e_sd, TRUE) +
                         log(pi) + fix_par[i, 'sigma0'])
              # accept or reject
              if(log(runif(1)) < r_d){
                phi[[i+1]] <- new_phi
                fix_par <- rbind(fix_par, new_fix)
              }
              else{
                phi[[i+1]] <- old_phi
                fix_par <- rbind(fix_par, fix_par[i,])
              }
            }
          }
          # move step
          else{
            old_phi <- phi[[i]]
            # choose a cluster to move
            idx <- sample(1:nrow(old_phi), size = 1)
            new_clust <- data.frame(x = rnorm(1, old_phi$x[idx], 0.01), y = rnorm(1, old_phi$y[idx], 0.01),
                                    N = rlnorm(1, log(old_phi$N[idx]), 0.05), R = rlnorm(1, log(old_phi$R[idx]), 0.05),
                                    e = rlnorm(1, log(old_phi$e[idx]), 0.05), n = rlnorm(1, log(old_phi$n[idx]), 0.05), 
                                    theta = rnorm(1, old_phi$theta[idx], 0.05), sigma = rlnorm(1, log(old_phi$sigma[idx]), 0.05))
            new_phi <- bind_rows(old_phi[-idx,], new_clust)
            # adaptive MCMC for Theta
            if(i <= 1000){
              new_b <- rnorm(1, fix_par[i, 'b0'], 0.025)
              old_b <- fix_par[i, 'b0']
              new_lc <- rnorm(1, fix_par[i, 'lc'], 0.025)
              old_lc <- fix_par[i, 'lc']
              new_s0 <- rnorm(1, fix_par[i, 'sigma0'], 0.01)
              old_s0 <- fix_par[i, 'sigma0']
              new_g <- rnorm(3*nrow(gal_start), fix_par[i, 4:(3+3*nrow(gal_start))], 0.025)
              old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
              new_fix <- c(new_b, new_lc, new_s0, new_g)
            }
            else if(i <= 10000){
              new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par))
              new_b <- new_fix[1]
              old_b <- fix_par[i, 'b0']
              new_lc <- new_fix[2]
              old_lc <- fix_par[i, 'lc']
              new_s0 <- new_fix[3]
              old_s0 <- fix_par[i, 'sigma0']
              new_g <- new_fix[4:(3+3*nrow(gal_start))]
              old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
            }
            else{
              new_fix <- MASS::mvrnorm(1, fix_par[i,], 0.1*cov(fix_par[1:10000,]))
              new_b <- new_fix[1]
              old_b <- fix_par[i, 'b0']
              new_lc <- new_fix[2]
              old_lc <- fix_par[i, 'lc']
              new_s0 <- new_fix[3]
              old_s0 <- fix_par[i, 'sigma0']
              new_g <- new_fix[4:(3+3*nrow(gal_start))]
              old_g <- fix_par[i, 4:(3+3*nrow(gal_start))]
            }
            
            new_gal <- matrix(new_g, nrow = nrow(gal_start), byrow = T)
            colnames(new_gal) <- c('N', 'R', 'n')
            old_gal <- matrix(old_g, nrow = nrow(gal_start), byrow = T)
            colnames(old_gal) <- c('N', 'R', 'n')
            
            # MH ratio
            r_m <- log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = new_gal, 
                                            par_UDG = new_phi, sigma0 = new_s0, b0 = new_b, 
                                            lc = new_lc, grid = grid, gal_prior = gal_prior, 
                                            UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                            lc_prior = lc_prior, s0_prior = s0_prior) -
                   log_post_superpose_color(X, Y, S, gal_fix = gal_fix, gal_rand = old_gal, 
                                            par_UDG = old_phi, sigma0 = old_s0, b0 = old_b, 
                                            lc = old_lc, grid = grid, gal_prior = gal_prior, 
                                            UDG_prior = UDG_prior, b0_prior = b0_prior, 
                                            lc_prior = lc_prior, s0_prior = s0_prior) +
                   dlnorm(old_phi$N[idx], log(new_clust$N), 0.05, TRUE) - dlnorm(new_clust$N, log(old_phi$N[idx]), 0.05, TRUE) +
                   dlnorm(old_phi$R[idx], log(new_clust$R), 0.05, TRUE) - dlnorm(new_clust$R, log(old_phi$R[idx]), 0.05, TRUE) +
                   dlnorm(old_phi$n[idx], log(new_clust$n), 0.05, TRUE) - dlnorm(new_clust$n, log(old_phi$n[idx]), 0.05, TRUE) +
                   dlnorm(old_phi$e[idx], log(new_clust$e), 0.05, TRUE) - dlnorm(new_clust$e, log(old_phi$e[idx]), 0.05, TRUE) +
                   dlnorm(old_phi$sigma[idx], log(new_clust$sigma), 0.05, TRUE) - dlnorm(new_clust$sigma, log(old_phi$sigma[idx]), 0.05, TRUE)
            # accept or reject
            if(log(runif(1)) < r_m){
              phi[[i+1]] <- new_phi
              fix_par <- rbind(fix_par, new_fix)
            }
            else{
              phi[[i+1]] <- old_phi
              fix_par <- rbind(fix_par, fix_par[i,])
            }
          }
        }
      }
      return(list(phi, fix_par))  
    }
  }
}
