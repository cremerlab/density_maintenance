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
    vector<lower=0>[N] volume = pi() .* widths^2 .* (3 .* lengths - widths);
    vector<lower=0>[N] SA = pi() .* widths .* lengths;
    vector<lower=0>[N] SAV = SA ./ volume;
}

parameters {
    // Hyper parameters
    real<lower=0> width_mu;
    real<lower=0> length_mu;
    real<lower=0> vol_mu;
    real<lower=0> SAV_mu;
    real tau;

    // Lower-level parameters
    vector[J] width_mu_1_tilde;
    vector[J] length_mu_1_tilde;
    vector<lower=0>[J] length_beta_1;
    vector[J] vol_mu_1_tilde;
    vector[J] sav_mu_1_tilde;;

    // Homoscedastic error parameters
    vector<lower=0>[J] homosced_width_sigma;
    vector<lower=0>[J] homosced_vol_sigma;
    vector<lower=0>[J] homosced_sav_sigma;
}

transformed parameters {
    // Perform uncentering
    vector<lower=0>[J] width_mu_1 = width_mu + tau * width_mu_1_tilde;
    vector<lower=0>[J] length_mu_1 = length_mu + tau * length_mu_1_tilde;
    vector<lower=0>[J] vol_mu_1 = vol_mu + tau * vol_mu_1_tilde;
    vector<lower=0>[J] sav_mu_1 = SAV_mu + tau * sav_mu_1_tilde;
    vector<lower=0>[J] length_alpha = length_beta_1 .* length_mu_1;
}


model {
    // Hyperparameters
    width_mu ~ std_normal(); 
    length_mu ~ std_normal();
    vol_mu ~ std_normal();
    SAV_mu ~ std_normal();
    tau ~ std_normal(); 

    // Low-level parameter
    width_mu_1_tilde  ~ std_normal(); 
    length_mu_1_tilde ~ std_normal(); 
    vol_mu_1_tilde ~ std_normal(); 
    sav_mu_1_tilde ~ std_normal(); 
    length_beta_1 ~ normal(0, 5);

    // Homoscedastic error parameters 
    homosced_width_sigma ~ std_normal();
    homosced_vol_sigma ~ std_normal();
    homosced_sav_sigma ~ std_normal();

    // Likelihoood
    widths ~ normal(width_mu_1[idx], homosced_width_sigma[idx]);
    lengths ~ gamma(length_alpha[idx], length_beta_1[idx]);   
    volume  ~ normal(vol_mu_1[idx], homosced_vol_sigma[idx]); 
    SAV  ~ normal(sav_mu_1[idx], homosced_sav_sigma[idx]); 
}
