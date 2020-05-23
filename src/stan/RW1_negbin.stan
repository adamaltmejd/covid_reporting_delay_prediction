data {
    int<lower=0> N;          //
    int Y[N];               //
   // real tau;
}
parameters {
    real<lower=0> tau;
    real<lower=0> phi;
    vector[N] theta;
}
model {
    tau ~ gamma(0.01,0.01);
    phi ~ gamma(0.001,0.001);

    for(i in 2:(N)){
    // non-centered parameterization
    //theta[i] ~ normal(2*theta[i-1]-theta[i-2], sqrt(1/tau));
    theta[i] ~ normal(theta[i-1], sqrt(1/tau));
    }
    Y ~ neg_binomial_2_log(theta, phi);
}
