data {
  int N;
  int d[N]; 
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
    int d_ppc[N];
    #vector[N] d_log_lik;
    for (n in 1:N) {
        d_ppc[n] = poisson_rng(rate);
        #d_log_lik[n] = poisson_lpmf(d[n]| rate);
    }
}