data {
    int<lower=1> N;
    vector<lower=0>[N] widths;
    vector<lower=0>[N] lengths;
}

parameters {
    real<lower=0> width_mu;
    real<lower=0> width_sigma;
    real<lower=0> length_mu;
    real<lower=0> length_sigma;
}

model {
    width_mu ~ normal(0, 1);
    width_sigma ~ normal(0, 0.1);
    length_mu ~ normal(0, 3);
    length_sigma ~ normal(0, 1);    

    widths ~ normal(width_mu, width_sigma);
    lengths ~ normal(length_mu, length_sigma);
}