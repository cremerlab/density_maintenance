data {
    // Dimensional information
    int<lower=1> J; // Number of biological replicates
    int<lower=1> N; // Number of cell measurements
    int<lower=1, upper=J> idx[N]; // ID vector mapping each cell to biological replicate
    real<lower=0> periplasmic_diam;

    // Measured information
    vector<lower=0>[N] widths;
    vector<lower=0>[N] lengths;
}

transformed data {
    // vector<lower=0>[N] aspect_ratio = lengths ./ widths;
    vector<lower=0>[N] volume = pi() .* ( widths .^ 2 ./ 12) .* (3 .* lengths - widths);
    // vector<lower=0>[N] surface_area = pi() .* lengths .* widths;
}

parameters {
    // Hyper parameters
    real<lower=0> width_mu;
    real<lower=0> length_mu;
    real<lower=0> vol_mu;
    // real<lower=0> sa_mu;
    // real<lower=0> ar_mu;
    real tau;

    // real<lower=0> width_sigma;
    // real<lower=0> length_sigma;
    // real<lower=0> vol_sigma;
    // real<lower=0> sa_sigma;
    // real<lower=0> ar_sigma;

    // Lower-level parameters
    vector[J] width_mu_1_tilde;
    vector[J] length_mu_1_tilde;
    vector[J] vol_mu_1_tilde;
    // vector[J] ar_mu_1_tilde;
    // vector[J] sa_mu_1_tilde;
    vector<lower=0>[J] homosced_width_sigma;
    vector<lower=0>[J] homosced_length_sigma;
    vector<lower=0>[J] homosced_vol_sigma;
    // vector<lower=0>[J] homosced_ar_sigma;
    // vector<lower=0>[J] homosced_sa_sigma;
}

transformed parameters {
    vector<lower=0>[J] width_mu_1 = width_mu + tau * width_mu_1_tilde;
    vector<lower=0>[J] length_mu_1 = length_mu + tau * length_mu_1_tilde;
    vector<lower=0>[J] vol_mu_1 = vol_mu + tau * vol_mu_1_tilde;
    // vector[J] ar_mu_1 = ar_mu + tau * ar_mu_1_tilde;
    // vector[J] sa_mu_1 = sa_mu + tau * sa_mu_1_tilde;
}

model {
    width_mu ~ gamma(4.2, 4.5);
    length_mu ~ gamma(4.5, 2.75);
    vol_mu ~ gamma(4.5, 2.75); 
    // sa_mu ~ std_normal();
    // ar_mu ~ normal(0, 2);

    // width_sigma ~ std_normal();
    // length_sigma ~ std_normal();
    // vol_sigma ~ std_normal();
    // sa_sigma ~ std_normal();
    // ar_sigma ~ std_normal();
          
    tau ~ normal(0, 1); 
    homosced_width_sigma ~ std_normal();
    homosced_length_sigma ~ std_normal();
    homosced_vol_sigma ~ std_normal();
    // homosced_sa_sigma ~ std_normal();
    // homosced_ar_sigma ~ std_normal();

    width_mu_1_tilde  ~ std_normal(); 
    length_mu_1_tilde ~ std_normal(); 
    vol_mu_1_tilde ~ std_normal(); 
    // sa_mu_1_tilde ~ std_normal(); 
    // ar_mu_1_tilde ~ std_normal();
    
    widths ~ normal(width_mu_1[idx], homosced_width_sigma[idx]);
    lengths ~ normal(length_mu_1[idx], homosced_length_sigma[idx]);   
    volume  ~ normal(vol_mu_1[idx], homosced_vol_sigma[idx]);
    // surface_area  ~ normal(sa_mu_1[idx], homosced_sa_sigma[idx]);
    // aspect_ratio  ~ normal(ar_mu_1[idx], homosced_ar_sigma[idx]);
 
}
