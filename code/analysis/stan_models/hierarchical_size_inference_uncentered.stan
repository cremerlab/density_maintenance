data {
    // Dimensional information
    int<lower=1> J; // Number of biological replicates
    int<lower=1> N; // Number of cell measurements
    array[N] int<lower=1, upper=J> idx; // ID vector mapping each cell to biological replicate

    // Measured information
    vector<lower=0>[N] widths;
    vector<lower=0>[N] lengths;
}

transformed data {
    vector<lower=0>[N] volume = pi() .* widths^2 .* (3 .* lengths - widths) / 12;
    vector<lower=0>[N] SA = pi() .* widths .* lengths;
    vector<lower=0>[N] SAV = SA ./ volume;
}

parameters {
    // Hyper parameters
    real<lower=0, upper=10> width_mu;
    real<lower=0> length_alpha;
    real<lower=0> length_beta;
    real<lower=0, upper=10> vol_mu;
    real<lower=0, upper=10> SAV_mu;
    real tau;

    // Lower-level parameters
    vector[J] width_mu_1_tilde;
    vector[J] length_beta_1_tilde;
    vector[J] length_alpha_1_tilde;
    vector[J] vol_mu_1_tilde;
    vector[J] sav_mu_1_tilde;;

    // Homoscedastic error parameters
    vector<lower=0>[J] homosced_width_sigma;
    vector<lower=0>[J] homosced_vol_sigma;
    vector<lower=0>[J] homosced_sav_sigma;
}

transformed parameters {
    // Perform uncentering
    vector<lower=0, upper=10>[J] width_mu_1 = width_mu + tau * width_mu_1_tilde;
    vector<lower=0, upper=10>[J] vol_mu_1 = vol_mu + tau * vol_mu_1_tilde;
    vector<lower=0, upper=10>[J] sav_mu_1 = SAV_mu + tau * sav_mu_1_tilde;
    vector<lower=0>[J] length_alpha_1 = length_alpha + tau * length_alpha_1_tilde;
    vector<lower=0>[J] length_beta_1 = length_beta + tau * length_beta_1_tilde ;
}


model {
    // Hyperparameters
    width_mu ~ gamma(6.5, 7.25); 
    length_alpha ~gamma(5.5, 2.25); 
    length_beta ~ gamma(5.5, 2.25);
    vol_mu ~ gamma(5.5, 2.25);
    SAV_mu ~ gamma(10, 2);
    tau ~ std_normal(); 

    // Low-level parameter
    width_mu_1_tilde  ~ std_normal(); 
    length_alpha_1_tilde ~ std_normal(); 
    length_beta_1_tilde ~ std_normal(); 
    vol_mu_1_tilde ~ normal(0, 3); 
    sav_mu_1_tilde ~ normal(0, 3); 

    // Homoscedastic error parameters 
    homosced_width_sigma ~ std_normal();
    homosced_vol_sigma ~ std_normal();
    homosced_sav_sigma ~ std_normal();

    // Likelihoood
    widths ~ normal(width_mu_1[idx], homosced_width_sigma[idx]);
    lengths ~ gamma(length_alpha_1[idx], length_beta_1[idx]);   
    volume  ~ normal(vol_mu_1[idx], homosced_vol_sigma[idx]); 
    SAV  ~ normal(sav_mu_1[idx], homosced_sav_sigma[idx]); 
}

generated quantities {
    vector[N] width_rep;
    vector[N] length_rep;
    vector[N] volume_rep;
    vector[N] SAV_rep;
    for (i in 1:N) {
        width_rep[i] = normal_rng(width_mu_1[idx[i]], homosced_width_sigma[idx[i]]);
        length_rep[i] = gamma_rng(length_alpha_1[idx[i]], length_beta_1[idx[i]]);
        volume_rep[i] = normal_rng(vol_mu_1[idx[i]], homosced_vol_sigma[idx[i]]);
        SAV_rep[i] = normal_rng(sav_mu_1[idx[i]], homosced_sav_sigma[idx[i]]);
    }
}
