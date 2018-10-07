data {
  int<lower=1> N;
  real<lower=-pi(),upper=pi()> rerr[N];
  int<lower=1,upper=2> cue[N];
  int<lower=1> J;
  int<lower=1,upper=J> id[N];
}

parameters {
  // group-level
  vector[2] mu;         // bias
  vector[2] logk;       // k (log scale)
  vector[2] logilambda;   // mixing parameter (probit scale)
  
  // subject specific
  vector<lower=0>[6] sigma_u;    // random effects standard deviations
  cholesky_factor_corr[6] L_u;   // L_u is the Choleski factor
  matrix[6,J] z_u;               // random effect matrix
}

transformed parameters {
  matrix<lower=0,upper=1>[2,J] lambda_j;
  matrix<lower=0>[2,J] kappa_j;
  matrix[2,J] mu_j;
  matrix[6,J] u;
  
  u = diag_pre_multiply(sigma_u, L_u) * z_u;
  
  for (j in 1:J){
    mu_j[1,j]= mu[1] + u[1,j];
    mu_j[2,j]= mu[2] + u[2,j];
    kappa_j[1,j]= exp(logk[1] + u[3,j]);
    kappa_j[2,j]= exp(logk[2] + u[4,j]);
    lambda_j[1,j]= inv_logit(logilambda[1] + u[5,j]);
    lambda_j[2,j]= inv_logit(logilambda[2] + u[6,j]);
  }
}

model {
  // priors
  L_u ~ lkj_corr_cholesky(2); 
  to_vector(z_u) ~ normal(0,1); 
  sigma_u ~ normal(0, 10);
  mu ~ normal(0.5, 10);
  logk ~ normal(2, 20);
  logilambda ~ normal(-4, 30);
 
  // likelihood
  for (n in 1:N){
    if (kappa_j[cue[n],id[n]] < 100)
      target += log_mix(lambda_j[cue[n],id[n]],
                      uniform_lpdf(rerr[n] | -pi(), pi()),
                      von_mises_lpdf(rerr[n] | mu_j[cue[n],id[n]], kappa_j[cue[n],id[n]]));
    else
      // for numerical stability, use normal approximation if vonMises's K > 100 (from Stan manual)
      target += log_mix(lambda_j[cue[n],id[n]],
                      uniform_lpdf(rerr[n] | -pi(), pi()),
                      normal_lpdf(rerr[n] | mu_j[cue[n],id[n]], sqrt(1/kappa_j[cue[n],id[n]])));
  }
}

