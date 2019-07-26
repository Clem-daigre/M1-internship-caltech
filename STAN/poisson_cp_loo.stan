data {
    int N;
    int d[N-1];
    int d_loo;  
    int n_loo; 
    real<lower=0> pes;
    real<lower=0> pls;
}

transformed data {
    real log_unif; 
    log_unif = -log(N-1);
}

parameters {
    real<lower=0> early;
    real<lower=0> late;
}

transformed parameters {
    vector[N-1] lp;
    {
        vector[N] lp_e;
        vector[N] lp_l;
        lp_e[1] = 0;
        lp_l[1] = 0;
        for (n in 1:(N-1)) {
            lp_e[n + 1] = lp_e[n] + poisson_lpmf(d[n] | early);
            lp_l[n + 1] = lp_l[n] + poisson_lpmf(d[n] | late);
        }
        lp = rep_vector(log_unif + lp_l[N], N-1) + head(lp_e, N-1) - head(lp_l, N-1);
    }
}


model {
    early ~ normal(0, pes);
    late ~ normal(0, pls);
    target += log_sum_exp(lp);
}


generated quantities {
    int<lower=1, upper=N-1> cp;
    real log_lik_d_loo;
    cp = categorical_rng(softmax(lp)); # effectively a gibbs sample of the changepoint? 
    if (cp < n_loo) {
        log_lik_d_loo = poisson_lpmf(d_loo | (n_loo < cp) ? early : late);
    }
    else {
        log_lik_d_loo = poisson_lpmf(d_loo | (n_loo < (cp+1)) ? early : late);
    }
}
