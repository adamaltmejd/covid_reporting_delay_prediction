data {
    int<lower=0> N;          //
    int Y[N];               //
   // real tau;
}
parameters {
    real<lower=0> tau;
    vector[N] theta;
}
model {
    tau ~ gamma(0.01,0.01);
    for(i in 2:(N)){
    // non-centered parameterization
    theta[i] ~ normal(theta[i-1], sqrt(1/tau));
    }
    Y ~ poisson_log(theta[1:N]);
}
