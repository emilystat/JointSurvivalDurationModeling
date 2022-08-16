library(mvtnorm)
library(rstan)
library(loo)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("function.R")

data1 <-
  read.csv("Emily_2009_data_updated_09062017.csv",
           header = T,
           na.strings = "NA")
N <- length(unique(data1[, 10]))
W <- matrix(0, nrow = N, ncol = N)
for (i in 1:N) {
  neigh.row <- which(data1[, 10] == i)[1]
  neigh <-
    as.numeric(unlist(strsplit(as.character(data1[neigh.row, 11]), ",")))
  W[i, neigh] <- 1
}
D <- diag(rowSums(W))

# set true model
model <- c("GMCAR", "UniCAR", "MvNorm", "IndNorm")
true.m <- model[3]

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
mcrun <- 50
niter <- 8000
nwarm <- 3000
nchains <- 4
setseed <- "2019"

eff.fitted1 <- matrix(nrow = mcrun, ncol = (11 * N + 6))
eff.fitted2 <- matrix(nrow = mcrun, ncol = (11 * N + 6))
eff.fitted3 <- matrix(nrow = mcrun, ncol = (11 * N + 6))
eff.fitted4 <- matrix(nrow = mcrun, ncol = (11 * N + 6))
eff.fitted5 <- matrix(nrow = mcrun, ncol = (11 * N + 6))
eff.true <- matrix(nrow = mcrun, ncol = (11 * N + 6))
ic.res <- matrix(nrow = mcrun, ncol = 35)
para.post <- matrix(nrow = mcrun, ncol = 21)

# parameters
rho <- 2
beta1 <- matrix(c(-12,-10,-8))
beta2 <- matrix(c(3, 4, 5))
alpha1 <- 0.7
alpha2 <- 0.9
tau1 <- 0.5
tau2 <- 2.5

# precision matrix
phi1sigma <- round(solve(tau1 * (D - alpha1 * W)), 6)
phi2sigma <- round(solve(tau2 * (D - alpha2 * W)), 6)

# more than one observations
group1 <- unlist(lapply(1:N, function(x)
  rep(x, 6)))
group2 <- unlist(lapply(1:N, function(x)
  rep(x, 3)))

mark1 <- Sys.time()
for (mc in 1:mcrun) {
  set.seed(2019 * mc)
  # random effect
  if (true.m == model[1]){
    # GMCAR
    eta0 = 0.5
    eta1 = 0.2
    phi1 <- t(rmvnorm(1, rep(0, N), phi1sigma, method = "chol"))
    phi2.mean <- (eta0 * diag(N) + eta1 * W) %*% (phi1)
    phi2 <- t(rmvnorm(1, phi2.mean, phi2sigma, method = "chol"))
  } else if (true.m == model[2]) {
    # seperate CAR
    phi1 <- t(rmvnorm(1, rep(0, N), phi1sigma, method = "chol"))
    phi2 <- t(rmvnorm(1, rep(0, N), phi2sigma, method = "chol"))
  } else if (true.m == model[3]) {
    # joint normal
    sigmamat <-
      rbind(cbind(2 * diag(N), 0.5 * diag(N)), cbind(0.5 * diag(N), 0.2 *
                                                       diag(N)))
    phi <- t(rmvnorm(1, c(rep(0, 2 * N)), sigmamat))
    phi1 <- phi[1:N, ]
    phi2 <- phi[-(1:N), ]
  } else if (true.m == model[4]) {
    # seperate normal
    sigma1 <- sqrt(2)
    sigma2 <- sqrt(0.2)
    phi1 <- as.matrix(rnorm(N, 0, sigma1))
    phi2 <- as.matrix(rnorm(N, 0, sigma2))
  }
  # summary(cbind(phi1, phi2))
  
  # fixed effect
  x2ind <- rep(1:3, N)
  x2 <- covar(x2ind)
  fixeff2 <- x2 %*% beta2
  
  # zero-truncated poisson
  y2 <- c()
  z2 <- c()
  for (i in 1:length(group2)) {
    y2[i] <- 0
    z2[i] <- phi2[group2[i]] + fixeff2[i]
    while (y2[i] == 0) {
      y2[i] <- rpois(1, exp(z2[i]))
    }
  }
  
  # fixed effect
  x1ind <- rep(1:3, N * 2)
  x1 <- covar(x1ind)
  fixeff1 <- x1 %*% beta1
  
  # Censored Weibull
  z1 <- c()
  O_uncensored <- c()
  for (i in 1:length(group1)) {
    z1[i] <- phi1[group1[i]] + fixeff1[i]
    O_uncensored[i] <-
      round(rweibull(1, shape = rho, scale = exp(-z1[i] / rho)))
  }
  
  # get censored data
  n_uncensor <- round(0.9 * length(O_uncensored))
  censor <- O_uncensored[order(O_uncensored)[n_uncensor]]
  
  # create data for stan
  censorno <- which(O_uncensored >= censor)
  group_censored <- group1[censorno]
  censor_time <- rep(censor, length(censorno))
  x1censor <- x1[censorno, ]
  uncensorno <- which(O_uncensored < censor)
  group_uncensored <- group1[uncensorno]
  t_uncensored <- O_uncensored[uncensorno]
  x1uncensor <- x1[uncensorno, ]
  
  # record true values
  eff.true[mc, ] <- c(beta1, beta2, phi1, phi2, z1, z2)
  
  # fit stan model
  full_d <-
    list(
      n = N,
      n_uncensor = length(uncensorno),
      n_censor = length(censorno),
      group_uncensored = group_uncensored,
      group_censored = group_censored,
      censor_time = censor_time,
      t_uncensored = t_uncensored,
      x1censor = x1censor,
      x1uncensor = x1uncensor,
      x2 = x2,
      n2 = length(group2),
      lower_limit = 1,
      group2 = group2,
      y2 = y2,
      W = W
    )
  
  full_fit1 <- stan(
    'GMCAR.stan',
    data = full_d,
    warmup = nwarm,
    iter = niter,
    chains = nchains,
    seed = setseed,
    verbose = FALSE,
    control = list(adapt_delta = 0.9)
  )
  eff.fitted1[mc, ] <- post.mean(full_fit1)
  ic.res[mc, 1:6] <- ic.fn(full_fit1)
  ic.res[mc, 7] <- div.fn(full_fit1, nchains = nchains)
  para.post[mc, 1:3] <- summary(full_fit1)$summary["eta0", c(4, 6, 8)]
  para.post[mc, 4:6] <- summary(full_fit1)$summary["eta1", c(4, 6, 8)]
  para.post[mc, 7:9] <- summary(full_fit1)$summary["tau1", c(4, 6, 8)]
  para.post[mc, 10:12] <- summary(full_fit1)$summary["tau2", c(4, 6, 8)]
  para.post[mc, 13:15] <- summary(full_fit1)$summary["alpha1", c(4, 6, 8)]
  para.post[mc, 16:18] <- summary(full_fit1)$summary["alpha2", c(4, 6, 8)]
  para.post[mc, 19:21] <- summary(full_fit1)$summary["rho", c(4, 6, 8)]
  
  full_fit2 <- stan(
    'GMCAR_rev.stan',
    data = full_d,
    warmup = nwarm,
    iter = niter,
    chains = nchains,
    seed = setseed,
    verbose = FALSE
  )
  eff.fitted2[mc, ] <- post.mean(full_fit2)
  ic.res[mc, 8:13] <- ic.fn(full_fit2)
  ic.res[mc, 14] <- div.fn(full_fit2, nchains = nchains)
  
  full_fit3 <- stan(
    'GMCAR_sep.stan',
    data = full_d,
    warmup = nwarm,
    iter = niter,
    chains = nchains,
    seed = setseed,
    verbose = FALSE
  )
  eff.fitted3[mc, ] <- post.mean(full_fit3)
  ic.res[mc, 15:20] <- ic.fn(full_fit3)
  ic.res[mc, 21] <- div.fn(full_fit3, nchains = nchains)
  
  full_fit4 <-
    stan(
      'Norm_joint.stan',
      data = full_d,
      warmup = nwarm,
      iter = niter,
      chains = nchains,
      seed = setseed,
      verbose = FALSE
    )
  eff.fitted4[mc, ] <- post.mean(full_fit4)
  ic.res[mc, 22:27] <- ic.fn(full_fit4)
  ic.res[mc, 28] <- div.fn(full_fit4, nchains = nchains)
  
  full_fit5 <-
    stan(
      'Norm_seperate.stan',
      data = full_d,
      warmup = nwarm,
      iter = niter,
      chains = nchains,
      seed = setseed,
      verbose = FALSE
    )
  eff.fitted5[mc, ] <- post.mean(full_fit5)
  ic.res[mc, 29:34] <- ic.fn(full_fit5)
  ic.res[mc, 35] <- div.fn(full_fit5, nchains = nchains)
  
  print(mc)
}

mark2 <- Sys.time()

amse <- matrix(nrow = 5, ncol = 6)
amse[1, ] <- amse.fun(eff.fitted1, eff.true)
amse[2, ] <- amse.fun(eff.fitted2, eff.true)
amse[3, ] <- amse.fun(eff.fitted3, eff.true)
amse[4, ] <- amse.fun(eff.fitted4, eff.true)
amse[5, ] <- amse.fun(eff.fitted5, eff.true)
colnames(amse) <- c("beta", "beta_se", "phi", "phi_se", "z", "z_se")

eff.fitted <- rbind(eff.fitted1, eff.fitted2, eff.fitted3, eff.fitted4, eff.fitted5)
eff.fitted <- cbind(eff.fitted, c(
  rep(1, mcrun),
  rep(2, mcrun),
  rep(3, mcrun),
  rep(4, mcrun),
  rep(5, mcrun)
))

if (true.m == model[1]){
  # GMCAR
  write.csv(eff.fitted, "GMCARhat.csv", row.names = FALSE)
  write.csv(eff.true, "GMCARtrue.csv", row.names = FALSE)
  write.csv(amse, "GMCARMSE.csv", row.names = c("GMCAR","reversed","UniCAR","MvNorm","IndNorm"))
  write.csv(ic.res, "GMCARIC.csv", row.names = FALSE)
  write.csv(para.post, "GMCARpara.csv", row.names = FALSE)
} else if (true.m == model[2]) {
  # seperate CAR
  write.csv(eff.fitted,"UniCARhat.csv",row.names = FALSE)
  write.csv(eff.true, "UniCARtrue.csv", row.names = FALSE)
  write.csv(amse,"UniCARMSE.csv",row.names = c("GMCAR","reversed","UniCAR","MvNorm","IndNorm"))
  write.csv(ic.res,"UniCARIC.csv",row.names = FALSE)
  write.csv(para.post, "UniCARpara.csv", row.names = FALSE)
} else if (true.m == model[3]) {
  # joint normal
  write.csv(eff.fitted,"MvNormhat.csv",row.names = FALSE)
  write.csv(eff.true, "MvNormtrue.csv", row.names = FALSE)
  write.csv(amse,"MvNormMSE.csv",row.names = c("GMCAR","reversed","UniCAR","MvNorm","IndNorm"))
  write.csv(ic.res,"MvNormIC.csv",row.names = FALSE)
  write.csv(para.post, "MvNormpara.csv", row.names = FALSE)
} else if (true.m == model[4]) {
  # seperate normal
  write.csv(eff.fitted,"IndNormhat.csv",row.names = FALSE)
  write.csv(eff.true, "IndNormtrue.csv", row.names = FALSE)
  write.csv(amse,"IndNormMSE.csv",row.names = c("GMCAR","reversed","UniCAR","MvNorm","IndNorm"))
  write.csv(ic.res,"IndNormIC.csv",row.names = FALSE)
  write.csv(para.post, "IndNormpara.csv", row.names = FALSE)
}

mark3 <- Sys.time()
print(mark3 - mark1)
