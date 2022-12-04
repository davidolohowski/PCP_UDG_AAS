source('src/help_functions.R')

BDM_MH <- function(X, S, beta_start, a_start, h_start, l_c_start, dx1, dx2, # likelihood parameters
                   bprior, aprior, hprior, lprior, M){
  
  phi <- list()
  fix_par <- matrix(c(beta_start, l_c_start), 1, 2)
  A <- (S[1,2] - S[1,1])*(S[2,2] - S[2,1])
  nc_start <- rpois(1, l_c_start*A) 
  
  phi[[1]] <- cbind(runif(nc_start, S[1,1], S[1,2]), runif(nc_start, S[2,1], S[2,2]), rep(a_start, nc_start), rep(h_start, nc_start))
  
  for (i in 1:M) {
    if(length(phi[[i]]) == 0){
      new_phi <- cbind(runif(1, S[1,1], S[1,2]), runif(1, S[2,1], S[2,2]), rlnorm(1, log(a_start), 0.5), rlnorm(1, log(h_start), 0.5))
      new_C <- rbind(phi[[i]][,1:2], new_phi[,1:2])
      old_C <- phi[[i]][,1:2]
      new_a <- c(phi[[i]][,3], new_phi[,3])
      old_a <- phi[[i]][,3]
      new_h <- c(phi[[i]][,4], new_phi[,4])
      old_h <- phi[[i]][,4]
      new_b <- rlnorm(1, log(fix_par[i, 1]), 0.05)
      old_b <- fix_par[i, 1]
      lc_new <- rlnorm(1, log(fix_par[i, 2]), 0.05)
      lc_old <- fix_par[i, 2]
      r_b <- log_post_SNCP(X, new_C, S, new_b, new_a, new_h, lc_new, dx1, dx2, bprior, aprior, hprior, lprior) + log(A) - 
        log_post_SNCP(X, old_C, S, old_b, old_a, old_h, lc_old, dx1, dx2, bprior, aprior, hprior, lprior) +
        dlnorm(lc_old, log(lc_new), 0.05, TRUE) - dlnorm(lc_new, log(lc_old), 0.05, TRUE) + 
        dlnorm(old_b, log(new_b), 0.05, TRUE) - dlnorm(new_b, log(old_b), 0.05, TRUE)
      if(log(runif(1)) < r_b){
        phi[[i+1]] <- rbind(phi[[i]], new_phi)
        fix_par <- rbind(fix_par, c(new_b, lc_new))
      }
      else{
        phi[[i+1]] <- phi[[i]]
        fix_par <- rbind(fix_par, c(old_b, lc_old))
      }
    }
    else{
      U_m <- runif(1)
      # birth or death
      if(U_m > 0.5){
        #birth
        U_b <- runif(1)
        if(U_b < 0.5){
          new_phi <- cbind(runif(1, S[1,1], S[1,2]), runif(1, S[2,1], S[2,2]), rlnorm(1, log(a_start), 0.5), rlnorm(1, log(h_start), 0.5))
          new_C <- rbind(phi[[i]][,1:2], new_phi[,1:2])
          old_C <- matrix(phi[[i]][,1:2], ncol = 2)
          new_a <- c(phi[[i]][,3], new_phi[,3])
          old_a <- phi[[i]][,3]
          new_h <- c(phi[[i]][,4], new_phi[,4])
          old_h <- phi[[i]][,4]
          new_b <- rlnorm(1, log(fix_par[i, 1]), 0.05)
          old_b <- fix_par[i, 1]
          lc_new <- rlnorm(1, log(fix_par[i, 2]), 0.05)
          lc_old <- fix_par[i, 2]
          r_b <- log_post_SNCP(X, new_C, S, new_b, new_a, new_h, lc_new, dx1, dx2, bprior, aprior, hprior, lprior) + log(A) - 
            log_post_SNCP(X, old_C, S, old_b, old_a, old_h, lc_old, dx1, dx2, bprior, aprior, hprior, lprior) - log(nrow(phi[[i]]) + 1) +
            dlnorm(lc_old, log(lc_new), 0.05, TRUE) - dlnorm(lc_new, log(lc_old), 0.05, TRUE) + 
            dlnorm(old_b, log(new_b), 0.05, TRUE) - dlnorm(new_b, log(old_b), 0.05, TRUE)
          if(log(runif(1)) < r_b){
            phi[[i+1]] <- rbind(phi[[i]], new_phi)
            fix_par <- rbind(fix_par, c(new_b, lc_new))
          }
          else{
            phi[[i+1]] <- phi[[i]]
            fix_par <- rbind(fix_par, c(old_b, lc_old))
          }
        }
        #death
        else{
          idx <- sample(1:nrow(phi[[i]]), size = 1)
          new_C <- matrix(phi[[i]][-idx,1:2], ncol = 2)
          old_C <- matrix(phi[[i]][,1:2], ncol = 2)
          new_a <- phi[[i]][-idx,3]
          old_a <- phi[[i]][,3]
          new_h <- phi[[i]][-idx,4]
          old_h <- phi[[i]][,4]
          new_b <- rlnorm(1, log(fix_par[i, 1]), 0.05)
          old_b <- fix_par[i, 1]
          lc_new <- rlnorm(1, log(fix_par[i, 2]), 0.05)
          lc_old <- fix_par[i, 2]
          r_d <- -1*(log_post_SNCP(X, old_C, S, old_b, old_a, old_h, lc_old, dx1, dx2, bprior, aprior, hprior, lprior) + log(A) - 
                       log_post_SNCP(X, new_C, S, new_b, new_a, new_h, lc_new, dx1, dx2, bprior, aprior, hprior, lprior) - log(nrow(phi[[i]]) + 1) -
                       dlnorm(lc_old, log(lc_new), 0.05, TRUE) + dlnorm(lc_new, log(lc_old), 0.05, TRUE) - 
                       dlnorm(old_b, log(new_b), 0.05, TRUE) + dlnorm(new_b, log(old_b), 0.05, TRUE))
          if(log(runif(1)) < r_d){
            phi[[i+1]] <- matrix(phi[[i]][-idx,], ncol = 4)
            fix_par <- rbind(fix_par, c(new_b, lc_new))
          }
          else{
            phi[[i+1]] <- phi[[i]]
            fix_par <- rbind(fix_par, c(old_b, lc_old))
          }
        }
      }
      # move
      else{
        idx <- sample(1:nrow(phi[[i]]), size = 1)
        new_phi <- cbind(rnorm(1, phi[[i]][idx,1], 0.05), rnorm(1, phi[[i]][idx,2], 0.05), rlnorm(1, log(phi[[i]][idx,3]), 0.05), rlnorm(1, log(phi[[i]][idx,4]), 0.05))
        new_C <- rbind(phi[[i]][-idx,1:2], new_phi[,1:2])
        old_C <- matrix(phi[[i]][,1:2], ncol = 2)
        new_a <- c(phi[[i]][-idx,3], new_phi[,3])
        old_a <- phi[[i]][,3]
        new_h <- c(phi[[i]][-idx,4], new_phi[,4])
        old_h <- phi[[i]][,4]
        new_b <- rlnorm(1, log(fix_par[i, 1]), 0.05)
        old_b <- fix_par[i, 1]
        lc_new <- rlnorm(1, log(fix_par[i, 2]), 0.05)
        lc_old <- fix_par[i, 2]
        r_m <- log_post_SNCP(X, new_C, S, new_b, new_a, new_h, lc_new, dx1, dx2, bprior, aprior, hprior, lprior) - 
          log_post_SNCP(X, old_C, S, old_b, old_a, old_h, lc_old, dx1, dx2, bprior, aprior, hprior, lprior) + 
          dlnorm(lc_old, log(lc_new), 0.05, TRUE) - dlnorm(lc_new, log(lc_old), 0.05, TRUE) + 
          dlnorm(old_b, log(new_b), 0.05, TRUE) - dlnorm(new_b, log(old_b), 0.05, TRUE) +
          dlnorm(phi[[i]][idx,3], log(new_phi[,3]), 0.05, TRUE) - dlnorm(new_phi[,3], log(phi[[i]][idx,3]), 0.05, TRUE) +
          dlnorm(phi[[i]][idx,4], log(new_phi[,4]), 0.05, TRUE) - dlnorm(new_phi[,4], log(phi[[i]][idx,4]), 0.05, TRUE)
        if(log(runif(1)) < r_m){
          phi[[i+1]] <- rbind(phi[[i]][-idx,], new_phi)
          fix_par <- rbind(fix_par, c(new_b, lc_new))
        }
        else{
          phi[[i+1]] <- phi[[i]]
          fix_par <- rbind(fix_par, c(old_b, lc_old))
        }
      }
    }
  }
  return(list(phi, fix_par))
}




































