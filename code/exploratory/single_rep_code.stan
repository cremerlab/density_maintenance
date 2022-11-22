 
functions {
    real elldiv_lpdf(vector cell_lengths, real ell_div, int N) { 
        vector[N] prob;
        for (i in 1:N) {
            // Determine if value should be updated or not
            prob[i] = ((cell_lengths[i] < (ell_div/2)) || (cell_lengths[i] > ell_div)) ?  1E-3 : (ell_div / (cell_lengths[i]^2)); 
        }
        return sum(log(prob));
    }
}
data { 
    // Define observables
    int<lower=1> N; // Number of measurements
    vector<lower=0>[N] cell_widths;
    vector<lower=0>[N] cell_lengths;
}

parameters { 
    //Instantiate the parameters
    real<lower=0> mu_omega;
    real<lower=0> sigma_omega;
    real<lower=0> alpha_ell;
    real<lower=0> beta_ell; 
}

transformed parameters {
    real<lower=0> mu_ell = alpha_ell / beta_ell;
}

model {
    // Priors
    //mu_omega ~ gamma(6, 4.3);
    mu_omega ~ std_normal();
    sigma_omega ~ std_normal(); 
    alpha_ell ~ normal(0, 5);
    beta_ell ~ normal(0, 5);

    // Likelihoods
    cell_widths ~ normal(mu_omega, sigma_omega);
    cell_lengths ~ gamma(alpha_ell, beta_ell);
}

generated quantities {
    vector[N] omega_rep ;
    vector[N] ell_rep;
    for (i in 1:N) { 
        omega_rep[i] = normal_rng(mu_omega, sigma_omega);
        ell_rep[i] = gamma_rng(alpha_ell, beta_ell);
}
}
