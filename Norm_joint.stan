data {
  int<lower = 1> n;
  int<lower=0> n_uncensor;
  int<lower=0> n_censor;
  int<lower = 1> n2;
  int<lower=1,upper=n> group_uncensored[n_uncensor];
  int<lower=1,upper=n> group_censored[n_censor];
  int<lower=1,upper=n> group2[n2];
  real<lower=0> t_uncensored[n_uncensor];
  real<lower=0> censor_time[n_censor];
  int lower_limit; 
  int<lower = lower_limit> y2[n2];
  matrix[n_uncensor,3] x1uncensor;
  matrix[n_censor,3] x1censor;
  matrix[n2,3] x2;
}
transformed data{
  matrix<lower = 0>[n, n] I;
  I = diag_matrix(rep_vector(1, n));
}
parameters {
  vector[3] beta1;
  vector[3] beta2;
  real<lower=0> rho;
  cov_matrix[2] sigma;
  vector[2*n] phi0;
  real<lower=1> t2_censored[n_censor];
}
transformed parameters {
  vector[n] phi1;
  vector[n] phi2;
  matrix[2,2] L;
  matrix[2,n] phi_noncenter;
  L = cholesky_decompose(sigma);
  for (i in 1:n){
  phi_noncenter[,i]= L * phi0[(2*i-1):(2*i)];
  }
  phi1 = phi_noncenter[1,]';
  phi2 = phi_noncenter[2,]';
}
model {
  phi0 ~ normal(0, 1);
  beta1 ~ normal(0, 100);
  beta2 ~ normal(0, 100);
  rho ~ gamma(1, 0.001);
  for(i in 1:n2){
    y2[i] ~ poisson(exp(phi2[group2[i]] + x2[i,]*beta2)) T[lower_limit, ];
  }
  for (i in 1:n_uncensor) {
	t_uncensored[i] ~ weibull(rho, exp(- phi1[group_uncensored[i]] / rho - x1uncensor[i,] * beta1 / rho));
  }
  for (i in 1:n_censor) {
	t2_censored[i] ~ weibull(rho, exp(- phi1[group_censored[i]] / rho - x1censor[i,] * beta1 / rho) / censor_time[i]);
  }
}
generated quantities {
  real log_lik[n2+n_uncensor+n_censor];
  for(i in 1:n2){
    log_lik[i] = poisson_lpmf(y2[i] | exp(phi2[group2[i]] + x2[i,] * beta2)) - poisson_lccdf(0 | exp(phi2[group2[i]] + x2[i,] * beta2));
  }
  for (i in 1:n_uncensor) {
	log_lik[n2+i] = weibull_lpdf(t_uncensored[i] | rho, exp(- phi1[group_uncensored[i]] / rho - x1uncensor[i,] * beta1 / rho));
  }
  for (i in 1:n_censor) {
	log_lik[n2+n_uncensor+i] = weibull_lccdf(censor_time[i] | rho, exp(- phi1[group_censored[i]] / rho - x1censor[i,] * beta1 / rho));
  }
}
