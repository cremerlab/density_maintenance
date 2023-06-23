data { 
    // Dimensional information
    int<lower=1> N_ms; // Number of mass spectrometry measurements
    int<lower=1> N_phiRb;
    int<lower=1> N_size; // Number of size measurements
    int<lower=1> N_ppc;  // Number of elements desired in ppcs. 
    int<lower=1> N_biomass;
    int<lower=0> linear_alpha;

    // Observed data
    vector<lower=0>[N_ms] m_peri; // Mass of periplasmic protein per cell
    vector<lower=1>[N_size] alpha;  // Length-to-width aspect ratio
    vector<lower=0>[N_size] size_lam;
    vector<lower=1>[N_biomass] rho_biomass;
    vector<lower=0>[N_phiRb] widths;
    vector<lower=0>[N_phiRb] phiRb;
    vector<lower=0>[N_phiRb] phiRb_lam;

    // Constants for calculations
    real<lower=0> beta_rp; // Conversion factor from R/P to phiRb
    real<lower=0> delta; // Periplasmic width
    vector<lower=0>[N_ppc] phiRb_range;
    vector<lower=0>[N_ppc] lam_range;
}

parameters {
    real<lower=0> rho_biomass_mu;
    real<lower=0> rho_biomass_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    array[1-linear_alpha] real<lower=1> alpha_mu; 
    array[linear_alpha] real<lower=1> alpha_int;
    array[linear_alpha] real<lower=0> alpha_slope;
    real<lower=0> alpha_sigma;
    real<lower=0> zeta;
    real<lower=0> width_sigma;
}

transformed parameters { 
    vector[N_size] alpha_mu_;
    vector[N_phiRb] phiRb_alpha_mu_;
    for (i in 1:N_size) {
        if (linear_alpha) {
            alpha_mu_[i] = alpha_int[1] + alpha_slope[1] * size_lam[i];
        }
        else {
            alpha_mu_[i]  = alpha_mu[1];
        }
    }
    for (i in 1:N_phiRb) {
        if (linear_alpha) {
            phiRb_alpha_mu_[i] = alpha_int[1] + alpha_slope[1] * phiRb_lam[i];
        }
        else { 
            phiRb_alpha_mu_[i] = alpha_mu[1];
        }
    }
}

model {
    // Priors
    rho_biomass_mu ~ normal(300, 100);
    rho_biomass_sigma ~ normal(0, 10);
    alpha_sigma ~ std_normal();
    m_peri_mu ~ normal(0, 10);
    m_peri_sigma ~ std_normal();
    zeta ~ normal(0, 20);
    if (linear_alpha) {
        alpha_int ~ std_normal();
        alpha_slope ~ std_normal();
    }
    else {
        alpha_mu ~ std_normal();
    }

    // Likelihoods
    widths ~ normal(zeta .* (1 + phiRb./beta_rp) .* 12 .* phiRb_alpha_mu_ ./ (rho_biomass_mu .* (3 .* phiRb_alpha_mu_ - 1)), width_sigma);
    alpha ~ normal(alpha_mu_, alpha_sigma); 
    rho_biomass ~ normal(rho_biomass_mu, rho_biomass_sigma);
    m_peri ~ normal(m_peri_mu, m_peri_sigma);

}


generated quantities {
    // PPcs
    vector<lower=0>[N_ppc] w_rep;
    vector<lower=0>[N_ppc] ell_rep;
    vector<lower=0>[N_ppc] vol_rep;
    vector<lower=0>[N_ppc] rho_peri_rep;
    vector<lower=0>[N_ppc] alpha_rep;
    vector<lower=0>[N_ppc] alpha_mu_0;
    vector<lower=0>[N_ppc] rho_biomass_rep;

    for (i in 1:N_ppc) {
        if (linear_alpha) { 
            alpha_rep[i] = normal_rng(alpha_int[1] + alpha_slope[1] * lam_range[i], alpha_sigma);
            alpha_mu_0[i] = alpha_int[1] + alpha_slope[1] * lam_range[i];
            }
        else {
            alpha_rep[i] = normal_rng(alpha_mu[1], alpha_sigma);
            alpha_mu_0[i] = alpha_mu[1];
            }
        rho_biomass_rep[i] = normal_rng(rho_biomass_mu, rho_biomass_sigma);
        w_rep[i] =  normal_rng(zeta * (1 + phiRb_range[i] ./ beta_rp) * 12 * alpha_mu_0[i] / (rho_biomass_mu * (3 * alpha_mu_0[i] - 1)), width_sigma);
        ell_rep[i] = alpha_mu_0[i] * w_rep[i];
        vol_rep[i] = pi() * w_rep[i]^3 * (3 * alpha_mu_0[i] - 1) / 12;
        rho_peri_rep[i] = m_peri_mu / (delta * pi() * alpha_mu_0[i] * w_rep[i]^2);
    }

}