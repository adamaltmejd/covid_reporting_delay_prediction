data {
    int<lower=0> N;          //
    int<lower=0> P;          //
    int Y[N];               //
    int n[N];               //
    matrix[N, P] X;
    real<lower =0> sigma_beta;
    real<lower =0> sigma_alpha;
}
parameters {
    vector[P] beta_alpha;
    vector[P] beta_beta;
}
model {
        vector[N] alpha = exp(X*beta_alpha);
        vector[N] beta = exp(X*beta_beta);
        beta_alpha  ~ normal(0, sigma_alpha);
        beta_beta   ~ normal(0, sigma_beta);

    Y ~ beta_binomial(n, alpha, beta);
}
