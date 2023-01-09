data {
    // Dimensional information
    int<lower=1> J; // Number of biological replicates
    int<lower=1> N; // Number of cell measurements
    array[N] int<lower=1, upper=J> idx; // ID vector mapping each cell to biological replicate
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
    // Hyper parameters for means
    real<lower=0> width_mu;
    real<lower=0> length_mu;
    real<lower=0> vol_mu;
    real<lower=0> SAV_mu;

    // Hyper parameters for sigmas
    real<lower=0> width_sigma;
    real<lower=0> length_sigma;
    real<lower=0> vol_sigma;
    real<lower=0> SAV_sigma;

    // Lower-level parameters
    vector<lower=0>[J] width_mu_1;
    vector<lower=0>[J] length_mu_1;
    vector<lower=0>[J] length_beta;
    vector<lower=0>[J] vol_mu_1;
    vector<lower=0>[J] sav_mu_1;;

    // Homoscedastic error parameters
    vector<lower=0>[J] homosced_width_sigma;
    vector<lower=0>[J] homosced_vol_sigma;
    vector<lower=0>[J] homosced_sav_sigma;
}

transformed parameters {
    vector<lower=0>[J] length_alpha = length_mu_1 .* length_beta;

}

model {
    // Hyperparameters
    width_mu ~ gamma(5, 5.5);
    length_mu ~ gamma(4.5, 4.5);
    vol_mu ~ gamma(4.5, 4.5);   
    SAV_mu ~ gamma(5, 1);
    width_sigma ~ std_normal();
    length_sigma ~ std_normal();
    vol_sigma ~ std_normal();
    SAV_sigma ~ std_normal();

    // Low-level parameters
    width_mu_1  ~ normal(width_mu, width_sigma);
    length_mu_1 ~ normal(length_mu, length_sigma); 
    length_beta ~ std_normal();
    vol_mu_1 ~ normal(vol_mu, vol_sigma); 
    sav_mu_1 ~ normal(SAV_mu, SAV_sigma); 

    // Homoscedastic error parameters 
    homosced_width_sigma ~ std_normal();
    homosced_vol_sigma ~ std_normal();
    homosced_sav_sigma ~ std_normal();

    // Likelihoood
    widths ~ normal(width_mu_1[idx], homosced_width_sigma[idx]);
    lengths ~ gamma(length_alpha[idx], length_beta[idx]);   
    volume  ~ normal(vol_mu_1[idx], homosced_vol_sigma[idx]); 
    SAV  ~ normal(sav_mu_1[idx], homosced_sav_sigma[idx]); 
}
