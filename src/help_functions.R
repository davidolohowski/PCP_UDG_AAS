UDG_ints_norm <- function(x, c, a, h){
  a*dnorm(x[,1], c[1], h)*dnorm(x[,2], c[2], h)
}

integrate_UDG_ints <- function(S, dx1, dx2, c, a, h){
  x1 <- seq(from = S[1,1], to = S[1,2], by = dx1)
  x2 <- seq(from = S[2,1], to = S[2,2], by = dx2)
  s <- expand.grid(x1, x2)
  return(sum(UDG_ints_norm(s, c, a, h))*dx1*dx2)
}

log_lik_SNCP <- function(X, C, S, beta, a, h, dx1, dx2){
  n_c <- nrow(C)
  A <- (S[1,2] - S[1,1])*(S[2,2] - S[2,1])
  loc_prod <- 0
  norm_const <- 0
  if(length(C) == 0){
    return(-beta*A + log(beta)*nrow(X))
  }
  else{
    for (i in 1:n_c) {
      norm_const <- norm_const + integrate_UDG_ints(S, dx1, dx2, C[i,], a[i], h[i])
      loc_prod <- loc_prod + UDG_ints_norm(X, C[i,], a[i], h[i])
    }
    return(-norm_const - beta*A + sum(log(beta + loc_prod)))
  }
}

log_lik_HPP <- function(C, l_c, S){
  n_c <- nrow(C)
  A <- (S[1,2] - S[1,1])*(S[2,2] - S[2,1])
  if(length(C) == 0){
    return(-l_c*A)
  }
  else{
    return(-l_c*A + n_c*log(l_c))
  }
}

log_post_SNCP <- function(X, C, S, beta, a, h, l_c, dx1, dx2, # likelihood parameters
                          bprior, aprior, hprior, lprior){
  n_c <- nrow(C)
  log_post <- log_lik_SNCP(X, C, S, beta, a, h, dx1, dx2) + log_lik_HPP(C, l_c, S) + 
    dlnorm(l_c, lprior[1], lprior[2], TRUE) + dlnorm(beta, bprior[1], bprior[2], TRUE)
  if(length(C) == 0){
    return(log_post)
  }
  else{  
    log_post <- log_post + sum(dlnorm(a, aprior[1], aprior[2], TRUE)) + sum(dlnorm(h, hprior[1], hprior[2], TRUE))
    return(log_post)
  }
}


Sersic_ints <- function(X, #point pattern
                        c, # center
                        N, # average number of points
                        R, # scale
                        e = 1, # ellipticity
                        n = 0.5, theta = 0 # Sersic index and inclination angle
){
  c1 <- as.numeric(c[1])
  c2 <- as.numeric(c[2])
  r <- (((X[,1] - c1)*cos(theta) - (X[,2]-c2)*sin(theta))^2/R^2 + ((X[,1] - c1)*sin(theta) + (X[,2]-c2)*cos(theta))^2/R^2/e^2)^0.5
  return(N*exp(-r^(1/n))/n/gamma(2*n)/2/pi/R^2/e)
}


integrate_Sersic <- function(grid, c, N, R, e = 1, n = 0.5, theta = 0){
  dx1 <- grid[[2]][1]
  dx2 <- grid[[2]][2]
  
  return(unname(sum(Sersic_ints(grid[[1]], c, N, R, e, n, theta))*dx1*dx2))
}


log_lik_superpose <- function(X, S, gal_fix = NULL, gal_rand = NULL, par_UDG, b0, grid){
  if(class(S)[1] != 'SpatialPolygons' & class(S)[1] != 'SpatialPolygonsDataFrame'){
    stop('Observational window S must be of class SpatialPolygons or SpatialPolygonsDataFrame.')
  }
  else{
    A <- sf::st_area(sf::st_as_sf(S))
    # likelihood when there is no normal galaxies or UDGs
    if(is.null(gal_fix) & nrow(par_UDG) == 0){
      return(-b0*A + log(b0)*nrow(X))
    }
    # likelihood when there are normal galaxies but no UDGs
    else if(!is.null(gal_fix) & nrow(par_UDG) == 0){
      if(nrow(gal_fix) != nrow(gal_rand)){
        stop('Number of normal galaxies should be the same for gal_fix and gal_rand.')
      }
      else{
        loc_prod <- 0
        norm_const <- 0
        for (i in 1:nrow(gal_fix)) {
          gal_c <- gal_fix[i,c('x','y')]
          gal_theta <- gal_fix[i,'theta']
          gal_e <- gal_fix[i,'e']
          
          gal_N <- gal_rand[i,'N']
          gal_R <- gal_rand[i,'R']
          gal_n <- gal_rand[i,'n']
          
          norm_const <- norm_const + integrate_Sersic(grid, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
          loc_prod <- loc_prod + Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
        }
        return(-norm_const - b0*A + sum(log(b0 + loc_prod)))
      }
    }
    # likelihood when there is no normal galaxies but there are  UDGs
    else if(is.null(gal_fix) & nrow(par_UDG) != 0){
      loc_prod <- 0
      norm_const <- 0
      for (i in 1:nrow(par_UDG)) {
        UDG_c <- par_UDG[i, c('x','y')]
        UDG_N <- par_UDG[i, 'N']
        UDG_R <- par_UDG[i, 'R']
        UDG_e <- par_UDG[i, 'e']
        UDG_n <- par_UDG[i, 'n']
        UDG_theta <- par_UDG[i, 'theta']
        norm_const <- norm_const + integrate_Sersic(grid, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
        loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      return(-norm_const - b0*A + sum(log(b0 + loc_prod)))
    }
    # likelihood when there are both normal galaxies and UDGs
    else{
      loc_prod <- 0
      norm_const <- 0
      if(nrow(gal_fix) != nrow(gal_rand)){
        stop('Number of normal galaxies should be the same for gal_fix and gal_rand.')
      }
      else{
        for (i in 1:nrow(gal_fix)) {
          gal_c <- gal_fix[i,c('x','y')]
          gal_theta <- gal_fix[i,'theta']
          gal_e <- gal_fix[i,'e']
          
          gal_N <- gal_rand[i,'N']
          gal_R <- gal_rand[i,'R']
          gal_n <- gal_rand[i,'n']
          
          norm_const <- norm_const + integrate_Sersic(grid, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
          loc_prod <- loc_prod + Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
        }
        for (i in 1:nrow(par_UDG)) {
          UDG_c <- par_UDG[i, c('x','y')]
          UDG_N <- par_UDG[i, 'N']
          UDG_R <- par_UDG[i, 'R']
          UDG_e <- par_UDG[i, 'e']
          UDG_n <- par_UDG[i, 'n']
          UDG_theta <- par_UDG[i, 'theta']
          norm_const <- norm_const + integrate_Sersic(grid, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
          loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
        }
        return(-norm_const - b0*A + sum(log(b0 + loc_prod)))
      }
    }
  }
}

log_lik_superpose_color <- function(X, Y, S, gal_fix = NULL, gal_rand = NULL, par_UDG, sigma0, b0, grid){
  if(class(S)[1] != 'SpatialPolygons' & class(S)[1] != 'SpatialPolygonsDataFrame'){
    stop('Observational window S must be of class SpatialPolygons or SpatialPolygonsDataFrame.')
  }
  else{
    A <- sf::st_area(sf::st_as_sf(S))
    # likelihood when there is no normal galaxies or UDGs
    if(is.null(gal_fix) & nrow(par_UDG) == 0){
      return(-b0*A + log(b0)*nrow(X) + sum(dnorm(Y, 0, sigma0, TRUE)))
    }
    # likelihood when there are normal galaxies but no UDGs
    else if(!is.null(gal_fix) & nrow(par_UDG) == 0){
      if(nrow(gal_fix) != nrow(gal_rand)){
        stop('Number of normal galaxies should be the same for gal_fix and gal_rand.')
      }
      else{
        loc_prod <- 0
        norm_const <- 0
        for (i in 1:nrow(gal_fix)) {
          gal_c <- gal_fix[i,c('x','y')]
          gal_theta <- gal_fix[i,'theta']
          gal_e <- gal_fix[i,'e']
          
          gal_N <- gal_rand[i,'N']
          gal_R <- gal_rand[i,'R']
          gal_n <- gal_rand[i,'n']
          
          norm_const <- norm_const + integrate_Sersic(grid, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
          loc_prod <- loc_prod + Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
        }
        return(-norm_const - b0*A + sum(log(b0 + loc_prod)) + sum(dnorm(Y, 0, sigma0, TRUE)))
      }
    }
    # likelihood when there is no normal galaxies but there are  UDGs
    else if(is.null(gal_fix) & nrow(par_UDG) != 0){
      loc_prod <- 0
      norm_const <- 0
      Y_lik <- as.data.frame(matrix(0, nrow(X), nrow(par_UDG)+1))
      p <- as.data.frame(matrix(0, nrow(X), nrow(par_UDG)))
      for (i in 1:nrow(par_UDG)) {
        UDG_c <- par_UDG[i, c('x','y')]
        UDG_N <- par_UDG[i, 'N']
        UDG_R <- par_UDG[i, 'R']
        UDG_e <- par_UDG[i, 'e']
        UDG_n <- par_UDG[i, 'n']
        UDG_theta <- par_UDG[i, 'theta']
        norm_const <- norm_const + integrate_Sersic(grid, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
        loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      for (i in 1:nrow(par_UDG)) {
        UDG_c <- par_UDG[i, c('x','y')]
        UDG_N <- par_UDG[i, 'N']
        UDG_R <- par_UDG[i, 'R']
        UDG_e <- par_UDG[i, 'e']
        UDG_n <- par_UDG[i, 'n']
        UDG_theta <- par_UDG[i, 'theta']
        UDG_sigma <- par_UDG[i, 'sigma']
        p[,i] <- Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)/(b0+loc_prod)
        Y_lik[,i+1] <- p[,i]*dnorm(Y, 0, UDG_sigma)
      }
      Y_lik[,1] <- (1 - rowSums(p))*dnorm(Y, 0, sigma0)
      return(-norm_const - b0*A + sum(log(b0 + loc_prod)) + sum(log(rowSums(Y_lik))))
    }
    # likelihood when there are both normal galaxies and UDGs
    else{
      loc_prod <- 0
      norm_const <- 0
      Y_lik <- as.data.frame(matrix(0, nrow(X), nrow(par_UDG)+1))
      p <- as.data.frame(matrix(0, nrow(X), nrow(par_UDG)))
      if(nrow(gal_fix) != nrow(gal_rand)){
        stop('Number of normal galaxies should be the same for gal_fix and gal_rand.')
      }
      else{
        for (i in 1:nrow(gal_fix)) {
          gal_c <- gal_fix[i,c('x','y')]
          gal_theta <- gal_fix[i,'theta']
          gal_e <- gal_fix[i,'e']
          
          gal_N <- gal_rand[i,'N']
          gal_R <- gal_rand[i,'R']
          gal_n <- gal_rand[i,'n']
          
          norm_const <- norm_const + integrate_Sersic(grid, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
          loc_prod <- loc_prod + Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
        }
        for (i in 1:nrow(par_UDG)) {
          UDG_c <- par_UDG[i, c('x','y')]
          UDG_N <- par_UDG[i, 'N']
          UDG_R <- par_UDG[i, 'R']
          UDG_e <- par_UDG[i, 'e']
          UDG_n <- par_UDG[i, 'n']
          UDG_theta <- par_UDG[i, 'theta']
          UDG_sigma <- par_UDG[i, 'sigma']
          norm_const <- norm_const + integrate_Sersic(grid, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
          loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
        }
        for (i in 1:nrow(par_UDG)) {
          UDG_c <- par_UDG[i, c('x','y')]
          UDG_N <- par_UDG[i, 'N']
          UDG_R <- par_UDG[i, 'R']
          UDG_e <- par_UDG[i, 'e']
          UDG_n <- par_UDG[i, 'n']
          UDG_theta <- par_UDG[i, 'theta']
          UDG_sigma <- par_UDG[i, 'sigma']
          p[,i] <- Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)/(b0+loc_prod)
          Y_lik[,i+1] <- p[,i]*dnorm(Y, 0, UDG_sigma)
        }
        Y_lik[,1] <- (1 - rowSums(p))*dnorm(Y, 0, sigma0)
        return(-norm_const - b0*A + sum(log(b0 + loc_prod)) + sum(log(rowSums(Y_lik))))
      }
    }
  }
}


log_lik_PPP <- function(Nc, lc, A){
  if(Nc == 0){
    return(-lc*A)
  }
  else{
    return(-lc*A + Nc*log(lc))
  }
}

log_post_superpose <- function(X, S,  gal_fix = NULL, gal_rand = NULL, par_UDG, b0, lc, grid,  # likelihood parameters
                               gal_prior = NULL, UDG_prior, b0_prior, lc_prior){
  A <- sf::st_area(sf::st_as_sf(S))
  Nc <- nrow(par_UDG)
  log_post <- log_lik_superpose(X, S, gal_fix, exp(gal_rand), par_UDG, exp(b0), grid) + log_lik_PPP(Nc, exp(lc), A) + 
    dnorm(lc, lc_prior[1], lc_prior[2], TRUE) + dnorm(b0, b0_prior[1], b0_prior[2], TRUE) + b0 + lc
  
  if(is.null(gal_rand) & Nc == 0){
    return(log_post)
  }
  
  else if(!is.null(gal_rand) & Nc == 0){
    for (i in 1:nrow(gal_rand)) {
      
      gal_N <- gal_rand[i,'N']
      gal_R <- gal_rand[i,'R']
      gal_n <- gal_rand[i,'n']
      
      gal_N_prior <- as.numeric(gal_prior[i, c('N_mean', 'N_sd')])
      gal_R_prior <- as.numeric(gal_prior[i, c('R_mean', 'R_sd')])
      gal_n_prior <- as.numeric(gal_prior[i, c('n_mean', 'n_sd')])
      
      log_post <- log_post + dnorm(gal_N, gal_N_prior[1], gal_N_prior[2], TRUE) + 
        dnorm(gal_R, gal_R_prior[1], gal_R_prior[2], TRUE) +
        dnorm(gal_n, gal_n_prior[1], gal_n_prior[2], TRUE) + 
        gal_N + gal_R + gal_n
    }
    return(log_post)
  }
  
  else if(is.null(gal_rand) & Nc != 0){
    
    UDG_N_prior <- as.numeric(UDG_prior[c('N_mean', 'N_sd')])
    UDG_R_prior <- as.numeric(UDG_prior[c('R_mean', 'R_sd')])
    UDG_e_prior <- as.numeric(UDG_prior[c('e_mean', 'e_sd')])
    UDG_n_prior <- as.numeric(UDG_prior[c('n_mean', 'n_sd')])
    
    for (i in 1:nrow(par_UDG)) {
      
      UDG_N <- par_UDG[i, 'N']
      UDG_R <- par_UDG[i, 'R']
      UDG_e <- par_UDG[i, 'e']
      UDG_n <- par_UDG[i, 'n']
      UDG_theta <- par_UDG[i, 'theta']
      
      log_post <- log_post + dlnorm(UDG_N, UDG_N_prior[1], UDG_N_prior[2], TRUE) +
        dlnorm(UDG_R, UDG_R_prior[1], UDG_R_prior[2], TRUE) + 
        dlnorm(UDG_e, UDG_e_prior[1], UDG_e_prior[2], TRUE) +
        dlnorm(UDG_n, UDG_n_prior[1], UDG_n_prior[2], TRUE) + 
        dunif(UDG_theta, 0, pi, TRUE)
    }
    return(log_post)
  }
  
  else{
    for (i in 1:nrow(gal_rand)){
      
      gal_N <- gal_rand[i,'N']
      gal_R <- gal_rand[i,'R']
      gal_n <- gal_rand[i,'n']
      
      gal_N_prior <- as.numeric(gal_prior[i, c('N_mean', 'N_sd')])
      gal_R_prior <- as.numeric(gal_prior[i, c('R_mean', 'R_sd')])
      gal_n_prior <- as.numeric(gal_prior[i, c('n_mean', 'n_sd')])
      
      log_post <- log_post + dnorm(gal_N, gal_N_prior[1], gal_N_prior[2], TRUE) + 
        dnorm(gal_R, gal_R_prior[1], gal_R_prior[2], TRUE) +
        dnorm(gal_n, gal_n_prior[1], gal_n_prior[2], TRUE) + 
        gal_N + gal_R + gal_n
    }
    
    UDG_N_prior <- as.numeric(UDG_prior[c('N_mean', 'N_sd')])
    UDG_R_prior <- as.numeric(UDG_prior[c('R_mean', 'R_sd')])
    UDG_e_prior <- as.numeric(UDG_prior[c('e_mean', 'e_sd')])
    UDG_n_prior <- as.numeric(UDG_prior[c('n_mean', 'n_sd')])
    
    for (i in 1:nrow(par_UDG)) {
      
      UDG_N <- par_UDG[i, 'N']
      UDG_R <- par_UDG[i, 'R']
      UDG_e <- par_UDG[i, 'e']
      UDG_n <- par_UDG[i, 'n']
      UDG_theta <- par_UDG[i, 'theta']
      
      log_post <- log_post + dlnorm(UDG_N, UDG_N_prior[1], UDG_N_prior[2], TRUE) +
        dlnorm(UDG_R, UDG_R_prior[1], UDG_R_prior[2], TRUE) + 
        dlnorm(UDG_e, UDG_e_prior[1], UDG_e_prior[2], TRUE) +
        dlnorm(UDG_n, UDG_n_prior[1], UDG_n_prior[2], TRUE) + 
        dunif(UDG_theta, 0, pi, TRUE)
    }
    return(log_post)
  }
}


log_post_superpose_color <- function(X, Y, S,  gal_fix = NULL, gal_rand = NULL, par_UDG, sigma0, b0, lc, grid,  # likelihood parameters
                                     gal_prior = NULL, UDG_prior, b0_prior, lc_prior, s0_prior){
  A <- sf::st_area(sf::st_as_sf(S))
  Nc <- nrow(par_UDG)
  log_post <- log_lik_superpose_color(X, Y, S, gal_fix, exp(gal_rand), par_UDG, exp(sigma0), exp(b0), grid) + log_lik_PPP(Nc, exp(lc), A) + 
    dnorm(lc, lc_prior[1], lc_prior[2], TRUE) + dnorm(b0, b0_prior[1], b0_prior[2], TRUE) + VGAM::dlgamma(sigma0, 0, s0_prior[1], s0_prior[2], TRUE)+
    b0 + lc + sigma0
  
  if(is.null(gal_rand) & Nc == 0){
    return(log_post)
  }
  
  else if(!is.null(gal_rand) & Nc == 0){
    for (i in 1:nrow(gal_rand)) {
      
      gal_N <- gal_rand[i,'N']
      gal_R <- gal_rand[i,'R']
      gal_n <- gal_rand[i,'n']
      
      gal_N_prior <- as.numeric(gal_prior[i, c('N_mean', 'N_sd')])
      gal_R_prior <- as.numeric(gal_prior[i, c('R_mean', 'R_sd')])
      gal_n_prior <- as.numeric(gal_prior[i, c('n_mean', 'n_sd')])
      
      log_post <- log_post + dnorm(gal_N, gal_N_prior[1], gal_N_prior[2], TRUE) + 
        dnorm(gal_R, gal_R_prior[1], gal_R_prior[2], TRUE) +
        dnorm(gal_n, gal_n_prior[1], gal_n_prior[2], TRUE) + 
        gal_N + gal_R + gal_n
    }
    return(log_post)
  }
  
  else if(is.null(gal_rand) & Nc != 0){
    
    UDG_N_prior <- as.numeric(UDG_prior[c('N_mean', 'N_sd')])
    UDG_R_prior <- as.numeric(UDG_prior[c('R_mean', 'R_sd')])
    UDG_e_prior <- as.numeric(UDG_prior[c('e_mean', 'e_sd')])
    UDG_n_prior <- as.numeric(UDG_prior[c('n_mean', 'n_sd')])
    
    for (i in 1:nrow(par_UDG)) {
      
      UDG_N <- par_UDG[i, 'N']
      UDG_R <- par_UDG[i, 'R']
      UDG_e <- par_UDG[i, 'e']
      UDG_n <- par_UDG[i, 'n']
      UDG_theta <- par_UDG[i, 'theta']
      UDG_sigma <- par_UDG[i, 'sigma']
      
      log_post <- log_post + dlnorm(UDG_N, UDG_N_prior[1], UDG_N_prior[2], TRUE) +
        dlnorm(UDG_R, UDG_R_prior[1], UDG_R_prior[2], TRUE) + 
        dlnorm(UDG_e, UDG_e_prior[1], UDG_e_prior[2], TRUE) +
        dlnorm(UDG_n, UDG_n_prior[1], UDG_n_prior[2], TRUE) + 
        dunif(UDG_theta, 0, pi, TRUE) + dunif(UDG_sigma, 0, exp(sigma0), TRUE)
    }
    return(log_post)
  }
  
  else{
    for (i in 1:nrow(gal_rand)){
      
      gal_N <- gal_rand[i,'N']
      gal_R <- gal_rand[i,'R']
      gal_n <- gal_rand[i,'n']
      
      gal_N_prior <- as.numeric(gal_prior[i, c('N_mean', 'N_sd')])
      gal_R_prior <- as.numeric(gal_prior[i, c('R_mean', 'R_sd')])
      gal_n_prior <- as.numeric(gal_prior[i, c('n_mean', 'n_sd')])
      
      log_post <- log_post + dnorm(gal_N, gal_N_prior[1], gal_N_prior[2], TRUE) + 
        dnorm(gal_R, gal_R_prior[1], gal_R_prior[2], TRUE) +
        dnorm(gal_n, gal_n_prior[1], gal_n_prior[2], TRUE) + 
        gal_N + gal_R + gal_n
    }
    
    UDG_N_prior <- as.numeric(UDG_prior[c('N_mean', 'N_sd')])
    UDG_R_prior <- as.numeric(UDG_prior[c('R_mean', 'R_sd')])
    UDG_e_prior <- as.numeric(UDG_prior[c('e_mean', 'e_sd')])
    UDG_n_prior <- as.numeric(UDG_prior[c('n_mean', 'n_sd')])
    
    for (i in 1:nrow(par_UDG)) {
      
      UDG_N <- par_UDG[i, 'N']
      UDG_R <- par_UDG[i, 'R']
      UDG_e <- par_UDG[i, 'e']
      UDG_n <- par_UDG[i, 'n']
      UDG_theta <- par_UDG[i, 'theta']
      UDG_sigma <- par_UDG[i, 'sigma']
      
      log_post <- log_post + dlnorm(UDG_N, UDG_N_prior[1], UDG_N_prior[2], TRUE) +
        dlnorm(UDG_R, UDG_R_prior[1], UDG_R_prior[2], TRUE) + 
        dlnorm(UDG_e, UDG_e_prior[1], UDG_e_prior[2], TRUE) +
        dlnorm(UDG_n, UDG_n_prior[1], UDG_n_prior[2], TRUE) + 
        dunif(UDG_theta, 0, pi, TRUE) + dunif(UDG_sigma, 0, exp(sigma0), TRUE)
    }
    return(log_post)
  }
}


PZ_X <- function(X, UDG_par, gal_par = NULL, b0){
  if(is.null(gal_par)){
    if(nrow(UDG_par) == 0){
      return(rep(0, nrow(X)))
    }
    else{
      loc_prod <- 0
      for (i in 1:nrow(UDG_par)) {
        UDG_c <- UDG_par[i, c('x','y')]
        UDG_N <- UDG_par[i, 'N']
        UDG_R <- UDG_par[i, 'R']
        UDG_e <- UDG_par[i, 'e']
        UDG_n <- UDG_par[i, 'n']
        UDG_theta <- UDG_par[i, 'theta']
        
        loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      return(1 - b0/(b0 + loc_prod))
    }
  }
  else{
    if(nrow(UDG_par) == 0){
      return(rep(0, nrow(X)))
    }
    else{
      loc_prod <- 0
      gb_ints <- 0
      for (i in 1:nrow(gal_par)) {
        gal_c <- gal_par[i,c('x','y')]
        gal_theta <- gal_par[i,'theta']
        gal_e <- gal_par[i,'e']
        gal_N <- gal_par[i,'N']
        gal_R <- gal_par[i,'R']
        gal_n <- gal_par[i,'n']
        
        gb_ints <- gb_ints + Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
      }
      for (i in 1:nrow(UDG_par)) {
        UDG_c <- UDG_par[i, c('x','y')]
        UDG_N <- UDG_par[i, 'N']
        UDG_R <- UDG_par[i, 'R']
        UDG_e <- UDG_par[i, 'e']
        UDG_n <- UDG_par[i, 'n']
        UDG_theta <- UDG_par[i, 'theta']
        
        loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      return(1 - (b0 + gb_ints)/(b0 + gb_ints + loc_prod))
    }
  }
}

PZ_XY <- function(X, Y, UDG_par, gal_par = NULL, b0, sigma0){
  if(is.null(gal_par)){
    if(nrow(UDG_par) == 0){
      return(rep(0, nrow(X)))
    }
    else{
      loc_prod <- 0
      for (i in 1:nrow(UDG_par)) {
        UDG_c <- UDG_par[i, c('x','y')]
        UDG_N <- UDG_par[i, 'N']
        UDG_R <- UDG_par[i, 'R']
        UDG_e <- UDG_par[i, 'e']
        UDG_n <- UDG_par[i, 'n']
        UDG_theta <- UDG_par[i, 'theta']
        UDG_sigma <- UDG_par[i, 'sigma']
        
        loc_prod <- loc_prod + dnorm(Y, 0, UDG_sigma)*Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      return(1 - b0*dnorm(Y, 0, sigma0)/(b0*dnorm(Y, 0, sigma0) + loc_prod))
    }
  }
  else{
    if(nrow(UDG_par) == 0){
      return(rep(0, nrow(X)))
    }
    else{
      loc_prod <- 0
      gb_ints <- 0
      for (i in 1:nrow(gal_par)) {
        gal_c <- gal_par[i,c('x','y')]
        gal_theta <- gal_par[i,'theta']
        gal_e <- gal_par[i,'e']
        gal_N <- gal_par[i,'N']
        gal_R <- gal_par[i,'R']
        gal_n <- gal_par[i,'n']
        
        gb_ints <- gb_ints + dnorm(Y, 0, sigma0)*Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
      }
      for (i in 1:nrow(UDG_par)) {
        UDG_c <- UDG_par[i, c('x','y')]
        UDG_N <- UDG_par[i, 'N']
        UDG_R <- UDG_par[i, 'R']
        UDG_e <- UDG_par[i, 'e']
        UDG_n <- UDG_par[i, 'n']
        UDG_theta <- UDG_par[i, 'theta']
        UDG_sigma <- UDG_par[i, 'sigma']
        
        loc_prod <- loc_prod + dnorm(Y, 0, UDG_sigma)*Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      return(1 - (b0*dnorm(Y, 0, sigma0) + gb_ints)/(b0*dnorm(Y, 0, sigma0) + gb_ints + loc_prod))
    }
  }
}

post_intesity <- function(X, UDG_par, gal_par = NULL, b0){
  if(is.null(gal_par)){
    if(nrow(UDG_par) == 0){
      return(rep(b0, nrow(X)))
    }
    else{
      loc_prod <- 0
      for (i in 1:nrow(UDG_par)) {
        UDG_c <- UDG_par[i, c('x','y')]
        UDG_N <- UDG_par[i, 'N']
        UDG_R <- UDG_par[i, 'R']
        UDG_e <- UDG_par[i, 'e']
        UDG_n <- UDG_par[i, 'n']
        UDG_theta <- UDG_par[i, 'theta']
        
        loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      return(b0 + loc_prod)
    }
  }
  else{
    if(nrow(UDG_par) == 0){
      gb_ints <- 0
      for (i in 1:nrow(gal_par)) {
        gal_c <- gal_par[i,c('x','y')]
        gal_theta <- gal_par[i,'theta']
        gal_e <- gal_par[i,'e']
        gal_N <- gal_par[i,'N']
        gal_R <- gal_par[i,'R']
        gal_n <- gal_par[i,'n']
        
        gb_ints <- gb_ints + Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
      }
      return(b0 + gb_ints)
    }
    else{
      loc_prod <- 0
      gb_ints <- 0
      for (i in 1:nrow(gal_par)) {
        gal_c <- gal_par[i,c('x','y')]
        gal_theta <- gal_par[i,'theta']
        gal_e <- gal_par[i,'e']
        gal_N <- gal_par[i,'N']
        gal_R <- gal_par[i,'R']
        gal_n <- gal_par[i,'n']
        
        gb_ints <- gb_ints + Sersic_ints(X, gal_c, gal_N, gal_R, gal_e, gal_n, gal_theta)
      }
      for (i in 1:nrow(UDG_par)) {
        UDG_c <- UDG_par[i, c('x','y')]
        UDG_N <- UDG_par[i, 'N']
        UDG_R <- UDG_par[i, 'R']
        UDG_e <- UDG_par[i, 'e']
        UDG_n <- UDG_par[i, 'n']
        UDG_theta <- UDG_par[i, 'theta']
        
        loc_prod <- loc_prod + Sersic_ints(X, UDG_c, UDG_N, UDG_R, UDG_e, UDG_n, UDG_theta)
      }
      return(b0 + gb_ints + loc_prod)
    }
  }
}

simulate_HPP <- function(b0, S){
  A <- sf::st_area(sf::st_as_sf(S))
  n <- rpois(1, b0*A)
  X <- as.data.frame(sp::spsample(S, n, type = 'random'))[,1:2]
  return(X)
}

simulate_Sersic <- function(N, c, R, n, e, theta){
  i <- 1
  xy <- numeric(0)
  if(N > 0){
    while (i <= N) {
      x <- runif(1, -30*R + c[1], 30*R + c[1])
      y <- runif(1, -30*R + c[2], 30*R + c[2])
      z <- runif(1, 0, 1/n/gamma(2*n)/2/pi/R^2/e)
      
      r <- (((x - c[1])*cos(theta) - (y - c[2])*sin(theta))^2/R^2 + ((x - c[1])*sin(theta) + (y - c[2])*cos(theta))^2/R^2/e^2)^0.5
      if (Sersic_ints(cbind(x,y), c, N, R, e, n, theta)/N > z){
        xy <- rbind(xy, c(x,y))
        i <- i+1
      }
    }
    return(xy)
  }
  else{
    return(xy)
  }
}

detect_metric <- function(U, R, prob, C_r, C_rc){
  x <- U[,1]
  y <- U[,2]
  n <- length(x)
  U_region <- list()
  
  for (i in 1:n) {
    Ur <- circle.polygon(x[i], y[i], radius = R[i], sides = 100, by.length = F, poly.type = "cartesian")
    Ur <- Polygon(Ur)
    Ur <- SpatialPolygons(list(Polygons(list(Ur),'region')))
    U_region[[i]] <- st_as_sf(SpatialPolygonsDataFrame(Ur, data.frame(id = Ur@polygons[[1]]@ID, row.names = Ur@polygons[[1]]@ID)))
  }
  
  HPD_nc <- lapply(lapply(HPDregionplot(as.mcmc(C_r), vars = c('x', 'y'), n = 400, h = c(bw.nrd(C_r$x), bw.nrd(C_r$y)), prob = prob), function(x) cbind(x$x, x$y)), function(x) Polygon(x))
  HPD_nc <- lapply(HPD_nc, function(x) SpatialPolygons(list(Polygons(list(x), paste0('region', x@labpt[1])))))
  HPD_nc <- SpatialPolygons(lapply(HPD_nc, function(x){x@polygons[[1]]}))
  HPD_nc <- st_as_sf(HPD_nc)
  HPD_nc$prob <- 1-prob
  HPD_nc$Model <- 'Model 1'
  
  HPD_c <- lapply(lapply(HPDregionplot(as.mcmc(C_rc), vars = c('x', 'y'), n = 400, h = c(bw.nrd(C_rc$x), bw.nrd(C_rc$y)), prob = prob), function(x) cbind(x$x, x$y)), function(x) Polygon(x))
  HPD_c <- lapply(HPD_c, function(x) SpatialPolygons(list(Polygons(list(x), paste0('region', x@labpt[1])))))
  HPD_c <- SpatialPolygons(lapply(HPD_c, function(x){x@polygons[[1]]}))
  HPD_c <- st_as_sf(HPD_c)
  HPD_c$prob <- 1-prob
  HPD_c$Model <- 'Model 2'
  
  U_area <- 0
  Acc_nc <- 0
  Acc_c <- 0
  FP_nc <- 0
  FP_c <- 0
  
  for (i in 1:n) {
    Acc_nc <- Acc_nc + abs((length(st_intersects(U_region[[i]], HPD_nc)[[1]]) != 0) - 1)
    Acc_c <- Acc_c + abs((length(st_intersects(U_region[[i]], HPD_c)[[1]]) != 0) - 1)
    U_area <- U_area + R[[i]]^2*pi
    FP_nc <- FP_nc + length(st_intersects(U_region[[i]], HPD_nc)[[1]])
    FP_c <- FP_c + length(st_intersects(U_region[[i]], HPD_c)[[1]])
  }
  
  Prec_nc <- sum(st_area(HPD_nc))/U_area
  Prec_c <- sum(st_area(HPD_c))/U_area
  
  FP_nc <- nrow(HPD_nc) - FP_nc
  FP_c <- nrow(HPD_c) - FP_c
  
  return(list(acc = Acc_nc - Acc_c, FP = FP_nc - FP_c, prec = Prec_nc/Prec_c, HPD_nc = HPD_nc, HPD_c = HPD_c))
}





























