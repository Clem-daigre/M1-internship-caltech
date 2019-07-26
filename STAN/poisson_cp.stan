data {
    int N;
    int d[N]
    real<lower=0> pes;
    real<lower=0> pls;
}

transformed data {
    real log_unif; 
    log_unif = -log(N);
}

parameters {
    real<lower=0> early;
    real<lower=0> late;
}

transformed parameters {
    vector[N] lp;
    {
        vector[N + 1] lp_e;
        vector[N + 1] lp_l;
        lp_e[1] = 0;
        lp_l[1] = 0;
        for (n in 1:N) {
            lp_e[n + 1] = lp_e[n] + poisson_lpmf(d[n] | early);
            lp_l[n + 1] = lp_l[n] + poisson_lpmf(d[n] | late);
        }
        lp = rep_vector(log_unif + lp_l[N+1], N) + head(lp_e, N) - head(lp_l, N);
    }
}

model {
    early ~ normal(0, pes);
    late ~ normal(0, pls);
    target += log_sum_exp(lp);
}

