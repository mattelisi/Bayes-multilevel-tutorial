data {
  int<lower=1> N;                  // number of observations
  real RT[N];                      // reaction time (transformed in seconds)
  int<lower=0,upper=9> Days[N];    // predictor
  int<lower=1> J;                  // number of subjects
  int<lower=1,upper=J> Subject[N]; // subject id
}

parameters {
  vector[2] beta;               // fixed-effects parameters
  real<lower=0> sigma_e;        // residual std
  vector<lower=0>[2] sigma_u;   // random effects standard deviations
  cholesky_factor_corr[2] L_u;  // L_u is the Choleski factor of the correlation matrix
  matrix[2,J] z_u;              // random effect matrix
}

transformed parameters {
  matrix[2,J] u;
  u = diag_pre_multiply(sigma_u, L_u) * z_u; // use Cholesky to set correlation
}

model {
  real mu; // conditional mean of the dependent variable

  //priors
  L_u ~ lkj_corr_cholesky(2);   // LKJ prior for the correlation matrix
  to_vector(z_u) ~ normal(0,1); // before Cholesky, random effects are normal variates with SD=1
  sigma_u ~ cauchy(0, 0.5);     // SD of random effects (vectorized)
  sigma_e ~ cauchy(0, 0.5);     // prior for residual standard deviation
  beta[1] ~ normal(0.3, 1);     // prior for fixed-effect intercept
  beta[2] ~ normal(0.01, 0.5);  // prior for fixed-effect slope

  //likelihood
  for (i in 1:N){
    mu = beta[1] + u[1,Subject[i]] + (beta[2] + u[2,Subject[i]])*Days[i];
    RT[i] ~ normal(mu, sigma_e);
  }
}

generated quantities {
  // simulate model for posterior predictive check
  real y_rep[N];
  
  matrix[2,J] u_rep; // simulate 'random effects'
    for (j in 1:J) {
    u_rep[1:2,j] = multi_normal_rng(rep_vector(0,2), diag_matrix(sigma_u) * (L_u*L_u') * diag_matrix(sigma_u));
  }
  
  for (i in 1:N){
    y_rep[i] = normal_rng(beta[1] + u_rep[1,Subject[i]] + (beta[2] + u_rep[2,Subject[i]])*Days[i], sigma_e);
  }
}
