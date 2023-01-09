data {
    // Dimensional parameters
    int<lower=1> N;

    // Observed parameters
    vector<lower=0>[N] elapsed_time;
    vector<lower=0>[N] optical_density;
}

transformed data {
    vector[N] log_optical_density = log(optical_density);
}

parameters {
    real<lower=0> mu;
    real<lower=0> sigma;
    real<lower=0> od_init;
}

transformed parameters {  
    real log_od_init = log(od_init);
}


model { 

    // Priors
    mu ~ std_normal();
    sigma ~ normal(0, 0.1); 
    od_init ~ normal(0, 0.1);

    // Likelihood
    log_optical_density ~ normal(log_od_init + mu .* elapsed_time, sigma);

}