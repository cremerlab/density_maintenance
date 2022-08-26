data {
    // Dimensional information
    int<lower=1> J; // Number of biological replicates
    int<lower=1> N; // Number of cell measurements
    int<lower=1, upper=J> idx[N]; // ID vector mapping each cell to biological replicate

    // Measured information
    vector<lower=0>[N] widths;
    vector<lower=0>[N] lengths;
}

parameters {
    // Hyper parameters
    real<lower=0> width_mu_0;
    real<lower=0> length_mu_0;
    real<lower=0> width_sigma_0;
    real<lower=0> length_sigma_0;

    // Lower level parameters
    vector<lower=0>[J] width_mu;
    vector<lower=0>[J] length_mu;

    // Homoscedastic error
    vector<lower=0>[J] width_sigma;
    vector<lower=0>[J] length_sigma;
}

model {
    // Hyper Priors
    width_mu_0 ~ normal(0, 1);
    width_sigma_0 ~ normal(0, 1);
    length_mu_0 ~ normal(0, 1);
    length_sigma_0 ~ normal(0, 1);

    // Level-0 priors
    width_mu ~ normal(width_mu_0, width_sigma_0);
    length_mu ~ normal(length_mu_0, length_sigma_0);
    width_sigma ~ normal(0, 1);  
    length_sigma ~ normal(0, 1);  

    widths ~ normal(width_mu[idx], width_sigma[idx]);
    lengths ~ normal(length_mu[idx], length_sigma[idx]);   
}
