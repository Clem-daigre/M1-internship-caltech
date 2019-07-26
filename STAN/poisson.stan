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