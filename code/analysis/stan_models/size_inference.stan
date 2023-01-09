data {
    // Dimensional information
    int<lower=1> N; // Number of cell measurements
    real<lower=0> periplasmic_diam;

    // Measured information
    vector<lower=0>[N] widths;
    vector<lower=0>[N] lengths;
}

transformed data {
    vector<lower=0>[N] volume = pi() .* ( (widths - 2 * periplasmic_diam).^ 2 ./ 12) .* (3 .* (lengths - 2 * periplasmic_diam) - (widths - 2 * periplasmic_diam));
    vector<lower=0>[N] SA = pi() .* widths .* lengths;
    vector<lower=0>[N] SAV = SA ./ volume;
}

parameters {
    // Hyper parameters
    real<lower=0> width_mu;
    real<lower=0> length_mu;
    real<lower=0> vol_mu; 
    real<lower=0> SAV_mu; 
    real<lower=0> homosced_width_sigma;
    real<lower=0> homosced_length_sigma;
    real<lower=0> homosced_vol_sigma;
    real<lower=0> homosced_sav_sigma;
}

model {
    width_mu ~ gamma(5, 5.5);
    length_mu ~ gamma(4.5, 4.5);
    vol_mu ~ gamma(4.5, 4.5);   
    SAV_mu ~ gamma(5, 1);
    homosced_width_sigma ~ std_normal();
    homosced_length_sigma ~ std_normal();
    homosced_vol_sigma ~ std_normal();
    homosced_sav_sigma ~ std_normal();

    widths ~ normal(width_mu, homosced_width_sigma);
    lengths ~ normal(length_mu, homosced_length_sigma);   
    volume  ~ normal(vol_mu, homosced_vol_sigma);
    SAV  ~ normal(SAV_mu, homosced_sav_sigma);
}
