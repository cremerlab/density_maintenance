data {
    // Dimensional information
    int<lower=1> N; // Number of cell measurements
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
       
    real<lower=0> homosced_width_sigma;
    real<lower=0> homosced_length_sigma;
    real<lower=0> homosced_vol_sigma;
    real<lower=0> homosced_ar_sigma;
    real<lower=0> homosced_sa_sigma;
}

model {
    width_mu ~ std_normal();
    length_mu ~ std_normal();
    vol_mu ~ std_normal();
    sa_mu ~ std_normal();
    ar_mu ~ normal(0, 5);

    homosced_width_sigma ~ std_normal();
    homosced_length_sigma ~ std_normal();
    homosced_vol_sigma ~ std_normal();
    homosced_sa_sigma ~ std_normal();
    homosced_ar_sigma ~ std_normal();

    widths ~ normal(width_mu, homosced_width_sigma);
    lengths ~ normal(length_mu, homosced_length_sigma);   
    volume  ~ normal(vol_mu, homosced_vol_sigma);
    surface_area  ~ normal(sa_mu, homosced_sa_sigma);
    aspect_ratio  ~ normal(ar_mu, homosced_ar_sigma); 
}
