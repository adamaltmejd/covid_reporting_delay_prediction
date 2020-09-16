data {
    int<lower=0> N;          //
    int Y[N];               //
   // real tau;
}
parameters {
    real<lower=0> tau;
    //real<lower=0> tau_x;
    real<lower=0> phi;
    vector[N+5] theta;
    //vector[N+5] x;
    //real<upper=1> alpha;
}
model {
    tau ~ gamma(01,0.01);
    //tau_x ~ gamma(0.01,0.01);
    phi ~ gamma(0.001,0.001);
    //alpha ~ normal(0,0.5);
    for(i in 2:(N+5)){
    // non-centered parameterization
    theta[i] ~ normal(theta[i-1], sqrt(tau));
    //x[i]     ~ normal(x[i-1], tau_x);
    //theta[i] ~ normal( (1-alpha) * theta[i-1] + alpha*x[i-1] , tau);
    }
    Y ~ neg_binomial_2_log(theta[1:N], phi);
}
