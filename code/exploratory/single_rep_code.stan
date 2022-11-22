 
data { 
    // Define observables
    int<lower=1> N; // Number of measurements
    vector<lower=0>[N] cell_widths;
    vector<lower=0>[N] cell_lengths;
}

transformed data {
    vector<lower=0>[N] cell_volumes = (pi()/12) * cell_widths.^2 .* (3 .* cell_lengths - cell_widths);
}

parameters { 
    //Instantiate the parameters
    real<lower=0> mu_omega;
    real<lower=0> sigma_omega;
    real<lower=0> alpha_ell;
    real<lower=0> beta_ell; 
    real<lower=0> mu_vol;
    real<lower=0> sigma_vol;
}

transformed parameters {
    real<lower=0> mu_ell = alpha_ell / beta_ell;
}

model {
    // Priors
    mu_omega ~ std_normal();
    sigma_omega ~ std_normal(); 
    mu_vol ~ std_normal();
    sigma_vol ~ std_normal();
    alpha_ell ~ normal(0, 5);
    beta_ell ~ normal(0, 5);

    // Likelihoods
    cell_widths ~ normal(mu_omega, sigma_omega);
    cell_lengths ~ gamma(alpha_ell, beta_ell);
    cell_volumes ~ normal(mu_vol, sigma_vol);
}

generated quantities {
    vector[N] omega_rep ;
    vector[N] ell_rep;
    vector[N] vol_rep;
    for (i in 1:N) { 
        omega_rep[i] = normal_rng(mu_omega, sigma_omega);
        ell_rep[i] = gamma_rng(alpha_ell, beta_ell);
        vol_rep[i] = normal_rng(mu_vol, sigma_vol);
    }
}
