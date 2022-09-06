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
    vector<lower=0>[N] cyto_width = widths - 2 .* periplasmic_diam;
    vector<lower=0>[N] aspect_ratio = lengths ./ widths;
    vector<lower=0>[N] volume = pi() .* ( cyto_width .^ 2 ./ 12) .* (3 .* lengths - cyto_width);
    vector<lower=0>[N] surface_area = pi() .* lengths .* widths;
}

parameters {
    // Hyper parameters
    real<lower=0> width_mu;
    real<lower=0> length_mu;
    real<lower=0> vol_mu;
    real<lower=0> sa_mu;
    real<lower=0> ar_mu;

    real<lower=0> width_sigma;
    real<lower=0> length_sigma;
    real<lower=0> vol_sigma;
    real<lower=0> sa_sigma;
    real<lower=0> ar_sigma;

    // Lower-level parameters
    vector<lower=0>[J] width_mu_1;
    vector<lower=0>[J] length_mu_1;
    vector<lower=0>[J] vol_mu_1;
    vector<lower=0>[J] ar_mu_1;
    vector<lower=0>[J] sa_mu_1;
    vector<lower=0>[J] homosced_width_sigma;
    vector<lower=0>[J] homosced_length_sigma;
    vector<lower=0>[J] homosced_vol_sigma;
    vector<lower=0>[J] homosced_ar_sigma;
    vector<lower=0>[J] homosced_sa_sigma;
}

model {
    width_mu ~ std_normal();
    length_mu ~ std_normal();
    vol_mu ~ std_normal();
    sa_mu ~ std_normal();
    ar_mu ~ normal(0, 5);

    if ( J > 1 )  
        width_sigma ~ std_normal();
        length_sigma ~ std_normal();
        vol_sigma ~ std_normal();
        sa_sigma ~ std_normal();
        ar_sigma ~ std_normal();
      

    homosced_width_sigma ~ std_normal();
    homosced_length_sigma ~ std_normal();
    homosced_vol_sigma ~ std_normal();
    homosced_sa_sigma ~ std_normal();
    homosced_ar_sigma ~ std_normal();

    width_mu_1 ~ normal(width_mu, width_sigma);
    length_mu_1 ~ normal(length_mu, length_sigma);
    vol_mu_1 ~ normal(vol_mu, vol_sigma);
    sa_mu_1 ~ normal(sa_mu, sa_sigma);
    ar_mu_1 ~ normal(ar_mu, ar_sigma);
    
    widths ~ normal(width_mu_1[idx], homosced_width_sigma[idx]);
    lengths ~ normal(length_mu_1[idx], homosced_length_sigma[idx]);   
    volume  ~ normal(vol_mu_1[idx], homosced_vol_sigma[idx]);
    surface_area  ~ normal(sa_mu_1[idx], homosced_sa_sigma[idx]);
    aspect_ratio  ~ normal(ar_mu_1[idx], homosced_ar_sigma[idx]);
 
}
