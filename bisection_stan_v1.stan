data {
  int<lower=1> N;
  int<lower=0,upper=1> rr[N];
  real ds[N];
  int<lower=0,upper=1> cond[N];
  int<lower=1> J;
  int<lower=1,upper=J> id[N];
  
  int<lower=1> N_ppc;   // input data for ppc
  real ds_ppc[N_ppc];
  int<lower=0,upper=1> cond_ppc[N_ppc];
}

parameters {
  vector[4] beta;                    // fixed-effects parameters
  vector<lower=0,upper=1>[J] lambda; // lapse rate
  vector<lower=0>[2] beta_ab;        // population parameters of lapse rate
  vector<lower=0>[4] sigma_u;        // random effects standard deviations
  cholesky_factor_corr[4] L_u;       // L_u is the Choleski factor of the correlation matrix
  matrix[4,J] z_u;                   // random effect matrix
}

transformed parameters {
  matrix[4,J] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u; // use Cholesky to set correlation
}

model {
  real mu; // linear predictor

  //priors
  L_u ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0,1); // before Cholesky, random effects are normal variates with SD=1
  sigma_u ~ cauchy(0, 1);       // SD of random effects (vectorized)
  beta[1] ~ normal(0, 2);       // prior for intercept (cond 1)
  beta[2] ~ normal(0, 4);       // prior for slope (cond 1)
  beta[3] ~ normal(0, 1);       // prior for diff. in intercept (cond 1)
  beta[4] ~ normal(0, 2);       // prior for diff. in slope (cond 1)
  beta_ab[1] ~ gamma(4,1);      // hyper priors for lapse rate distribution
  beta_ab[2] ~ gamma(1,0.2);
  lambda ~ beta(beta_ab[1] ,beta_ab[2]);  //lapse rate

  //likelihood
  for (i in 1:N){
    mu = beta[1] + u[1,id[i]] + cond[i]*(beta[2]+u[2,id[i]]) 
        + (beta[3] + u[3,id[i]] + cond[i]*(beta[4] + u[4,id[i]] ))*ds[i];
    rr[i] ~ bernoulli((1-lambda[id[i]])*Phi(mu)+lambda[id[i]]/2);
  }
}

generated quantities {
  real temp_mu;
  real y_rep[N_ppc];
  vector[4] u_rep; // simulate 'random effects'
  real<lower=0,upper=1> lambda_rep;
  
  u_rep = multi_normal_cholesky_rng(rep_vector(0, 4), diag_matrix(sigma_u) * (L_u));
  lambda_rep = beta_rng(1 ,4);
  
  for (i in 1:N_ppc){
    temp_mu = beta[1] + u_rep[1] + cond_ppc[i]*(beta[2]+u_rep[2]) 
        + (beta[3] + u_rep[3] + cond_ppc[i]*(beta[4] + u_rep[4] ))*ds_ppc[i];
    y_rep[i] = bernoulli_rng((1-lambda_rep)*Phi(temp_mu)+lambda_rep/2);
  }
}