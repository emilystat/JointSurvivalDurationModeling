# AMSE
amse.fun <- function(fitted, true) {
  beta_mat <- fitted[, 1:6]
  phi_mat <- fitted[, 7:(7 + 2 * N - 1)]
  z_mat <- fitted[, (7 + 2 * N):ncol(fitted)]
  
  beta_true <- true[, 1:6]
  phi_true <- true[, 7:(7 + 2 * N - 1)]
  z_true <- true[, (7 + 2 * N):ncol(true)]
  
  amse <- c()
  amse[1] <- mean((beta_mat - beta_true) ^ 2)
  amse[2] <-
    mean(((beta_mat - beta_true) ^ 2 - amse[1]) ^ 2) / (nrow(beta_mat) * ncol(beta_mat) - 1)
  amse[3] <- mean((phi_mat - phi_true) ^ 2)
  amse[4] <-
    mean(((phi_mat - phi_true) ^ 2 - amse[3]) ^ 2) / (nrow(phi_mat) * ncol(phi_mat) - 1)
  amse[5] <- mean((z_mat - z_true) ^ 2)
  amse[6] <-
    mean(((z_mat - z_true) ^ 2 - amse[5]) ^ 2) / (nrow(z_mat) * ncol(z_mat) - 1)
  return(amse)
}

# transform covariate into binary matrix
covar <- function(vec) {
  xmat <- matrix(nrow = length(vec), ncol = 3)
  for (i in 1:length(vec)) {
    if (vec[i] == 1) {
      xmat[i,] = c(1, 0, 0)
    } else if (vec[i] == 2) {
      xmat[i,] = c(0, 1, 0)
    } else{
      xmat[i,] = c(0, 0, 1)
    }
  }
  return(xmat)
}

# get total divergences of a stanfit object
div.fn <- function(stanfit, nchains) {
  div <-
    sum(get_sampler_params(stanfit, inc_warmup = FALSE)[[1]][, "divergent__"])
  for (i in 2:nchains) {
    div <-
      div + sum(get_sampler_params(stanfit, inc_warmup = FALSE)[[i]][, "divergent__"])
  }
  return(div)
}

# return DIC and LOOIC statistics
ic.fn <- function(stanfit, nthin = 1) {
  stan.mat <- as.matrix(stanfit)
  
  # posterior mean of parameters
  para.est <-
    apply(stan.mat[seq(1, nrow(stan.mat), by = nthin),], 2, mean)
  fit.beta2 <- para.est[c("beta2[1]", "beta2[2]", "beta2[3]")]
  fit.beta1 <- para.est[c("beta1[1]", "beta1[2]", "beta1[3]")]
  
  poilik <- c()
  uncensor.weilik <- c()
  censor.weilik <- c()
  for (j in 1:length(y2)) {
    phi2col <- paste("phi2[", group2[j], "]", sep = "")
    poilik[j] <-
      log(dpois(y2[j], exp(para.est[phi2col] + x2[j,] %*% fit.beta2)) /
            (1 - dpois(0, exp(
              para.est[phi2col] + x2[j,] %*% fit.beta2
            ))))
  }
  for (j in 1:length(group_uncensored)) {
    uncensor.phi1col <-
      paste("phi1[", group_uncensored[j], "]", sep = "")
    uncensor.weilik[j] <-
      dweibull(
        t_uncensored[j],
        shape = para.est["rho"],
        scale = exp(-para.est[uncensor.phi1col] / para.est["rho"]
                    - x1uncensor[j,] %*% fit.beta1 / para.est["rho"]),
        log = TRUE
      )
  }
  for (j in 1:length(group_censored)) {
    censor.phi1col <- paste("phi1[", group_censored[j], "]", sep = "")
    censor.weilik[j] <-
      pweibull(
        censor_time[j],
        shape = para.est["rho"],
        scale = exp(-para.est[censor.phi1col] / para.est["rho"]
                    - x1censor[j,] %*% fit.beta1 / para.est["rho"]),
        lower.tail = FALSE,
        log.p = TRUE
      )
  }
  
  # LOOIC part
  log.lik.mat0 <-
    extract_log_lik(stanfit,
                    parameter_name = "log_lik",
                    merge_chains = TRUE)
  log.lik.mat <-
    log.lik.mat0[seq(1, nrow(log.lik.mat0), by = nthin),]
  loo.res <- loo(log.lik.mat)
  loo <- loo.res$estimates[3,]
  
  # DIC part
  Dbar <- -2 * mean(apply(log.lik.mat, 1, sum))
  Dhat <-
    -2 * (sum(poilik) + sum(uncensor.weilik) + sum(censor.weilik))
  pD <- Dbar - Dhat
  DIC <- 2 * Dbar - Dhat
  
  ic.res <- c(Dbar, Dhat, pD, DIC, loo)
  names(ic.res) <- c("Dbar", "Dhat", "pD", "DIC", "LOO", "LOOse")
  return(ic.res)
}

# get center value for maps
midfun <- function(vec) {
  min(vec) + diff(range(vec)) / 2
}

# MSE
post.mean <- function(stanfit) {
  fit_para <- extract(stanfit)
  
  phi1_est <- fit_para$phi1
  phi1_hat <- apply(phi1_est, 2, mean)
  
  beta1_est <- fit_para$beta1
  beta1_hat <- apply(beta1_est, 2, mean)
  fix1_hat <- x1 %*% beta1_hat
  z1_hat <- c()
  for (i in 1:length(group1)) {
    z1_hat[i] <- phi1_hat[group1[i]] + fix1_hat[i]
  }
  
  phi2_est <- fit_para$phi2
  phi2_hat <- apply(phi2_est, 2, mean)
  
  beta2_est <- fit_para$beta2
  beta2_hat <- apply(beta2_est, 2, mean)
  fix2_hat <- x2 %*% beta2_hat
  z2_hat <- c()
  for (i in 1:length(group2)) {
    z2_hat[i] <- phi2_hat[group2[i]] + fix2_hat[i]
  }
  
  return(c(beta1_hat, beta2_hat, phi1_hat, phi2_hat, z1_hat, z2_hat))
}