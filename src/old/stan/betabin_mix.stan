data {
    int<lower=0> N;          //
    int<lower=0> P_mu;          //
    int<lower=0> P_M;          //
    int ndays;

    int Y[N];               //
    int n[N];               //
    int day[N];               //
    matrix[N, P_mu] X_mu;
    matrix[N, P_M] X_M;
    real<lower =0> sigma_mu;
    real<lower =0> sigma_M;
}
parameters {
    vector[P_mu] beta_mu;
    vector[P_M] beta_M;
    simplex[2] theta;          // mixing proportions
    real<lower=0> mu_add;
    real<lower=0> M_sub;
}
model {
    vector[N] mu = inv_logit(X_mu*beta_mu);
    vector[N] M = exp(X_M*beta_M);
    vector[N] M0 = exp(X_M*beta_M - M_sub);
    vector[N] alpha = M .* mu;
    vector[N] mu0 = inv_logit(X_mu*beta_mu + mu_add);
    vector[N] alpha_0 = M0 .* mu0;
    vector[N] beta  = (1-mu) .* M;
    vector[N] beta0  = (1-mu0) .* M0;
    vector[2] lps[ndays];

    for(i in 1:ndays)
    {
        lps[i]=log(theta);
    }
    beta_mu ~ normal(0, sigma_mu);
    beta_M  ~ normal(0, sigma_M);
    mu_add  ~ gamma(1.1,0.001);
    M_sub   ~ gamma(1.1,0.001);
    for(i in 1:N) {
        if(day[i]==0){
             target += beta_binomial_lpmf(Y[i]|n[i],alpha[i],beta[i]);
        }else{
         lps[day[i]][1] += beta_binomial_lpmf(Y[i]|n[i],alpha_0[i],beta0[i]);
         lps[day[i]][2] += beta_binomial_lpmf(Y[i]|n[i],alpha[i],beta[i]);
        }
    }
    for(i in 1:ndays)
        target += log_sum_exp(lps[i]);

}
