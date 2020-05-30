data {
    int<lower=0> N;          //
    int<lower=0> P_mu;          //
    int<lower=0> P_M;          //

    int Y[N];               //
    int n[N];               //
    matrix[N, P_mu] X_mu;
    matrix[N, P_M] X_M;
    real<lower =0> sigma_mu;
    real<lower =0> sigma_M;
}
parameters {
    vector[P_mu] beta_mu;

    vector[P_M] beta_M;
}
model {
    vector[N] mu_ =X_mu*beta_mu;
    vector[N] M_ = X_M*beta_M;
    beta_mu ~ normal(0, sigma_mu);
    beta_M  ~ normal(0, sigma_M);
    for(i in 1:N) {
        real alpha = inv_logit(mu_[i]) * exp(M_[i]);
        real beta = (1-inv_logit(mu_[i])) * exp(M_[i]);
         target += beta_binomial_lpmf(Y[i]|n[i],alpha,beta);
    }

}
