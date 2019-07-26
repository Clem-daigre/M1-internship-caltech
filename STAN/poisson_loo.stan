data {
  int N;
  int d[N-1];
  int d_loo;  
  int n_loo; 
  real<lower=0> prs;
}

parameters {
    real<lower=0> rate;
}

model {
    rate ~ normal(0, prs);
    d ~ poisson(rate);
}


generated quantities {
    real log_lik_d_loo;
    log_lik_d_loo = poisson_lpmf(d_loo | rate);
}