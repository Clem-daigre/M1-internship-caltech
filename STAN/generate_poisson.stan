data {
  int<lower=1> N;
  real<lower=0> rate;
}

generated quantities {
  int d[N] = rep_array(0, N);
  for (n in 1:N) {
        d[n] = poisson_rng(rate);
  }
}
