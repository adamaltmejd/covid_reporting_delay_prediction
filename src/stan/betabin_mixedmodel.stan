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
    vector[ndays] mu_day;
    //vector[ndays] M_day;
    //real<lower =0> sigma_M_day;
    real<lower =0> sigma_mu_day;

    vector[P_M] beta_M;
}
model {
    vector[N] mu_ =X_mu*beta_mu;
    vector[N] M_ = X_M*beta_M;
    vector[P_M] expbeta_M = exp(beta_M);
    vector[2] lps[ndays];
    sigma_mu_day ~ gamma(1,0.01);
    //sigma_M_day  ~ gamma(1,0.01);
    mu_day  ~ normal(0, sigma_mu_day);
    //M_day   ~ normal(0, sigma_M_day);
    beta_mu ~ normal(0, sigma_mu);
    target += exponential_lpdf(expbeta_M|  1);
    for(i in 1:N) {

        if(day[i]==0){
            real beta = (1-inv_logit(mu_[i])) * exp(M_[i]);
            real alpha = inv_logit(mu_[i]) * exp(M_[i]);

             target += beta_binomial_lpmf(Y[i]|n[i],alpha,beta);
        }else{
            real alpha = inv_logit(mu_[i]   + mu_day[day[i]]) * exp(M_[i]  );
            real beta = (1-inv_logit(mu_[i] + mu_day[day[i]])) * exp(M_[i] ); // + M_day[day[i]]
            target += beta_binomial_lpmf(Y[i]|n[i],alpha,beta);
        }
    }

}
