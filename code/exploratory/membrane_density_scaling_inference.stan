data { 
    // Dimensional information
    int<lower=1> N_ms; // Number of mass spectrometry measurements
    int<lower=1> N_phiRb;
    int<lower=1> N_size; // Number of size measurements
    int<lower=1> N_ppc;  // Number of elements desired in ppcs. 
    int<lower=1> N_biomass;

    // Observed data
    vector<lower=0>[N_ms] m_peri; // Mass of periplasmic protein per cell
    vector<lower=1>[N_size] alpha;  // Length-to-width aspect ratio
    vector<lower=1>[N_biomass] rho_biomass;
    vector<lower=0>[N_phiRb] widths;
    vector<lower=0>[N_phiRb] phiRb;

    // Constants for calculations
    real<lower=0> beta_rp; // Conversion factor from R/P to phiRb
    real<lower=0> delta; // Periplasmic width
    vector<lower=0>[N_ppc] phiRb_range;
}

parameters {
    real<lower=0> rho_biomass_mu;
    real<lower=0> rho_biomass_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    real<lower=1> alpha_mu; 
    real<lower=0> alpha_sigma;
    real<lower=0> zeta;
    real<lower=0> width_sigma;
}

transformed parameters { 
    real<lower=0> kappa = 12 * alpha_mu  / (3 * alpha_mu - 1);
    vector[N_phiRb] width_mu = zeta * kappa .* (1 + phiRb ./ beta_rp)  ./ rho_biomass_mu;
}

model {
    // Priors
    rho_biomass_mu ~ normal(300, 100);
    rho_biomass_sigma ~ normal(0, 10);
    alpha_mu ~ normal(1, 2);
    alpha_sigma ~ std_normal();
    m_peri_mu ~ normal(0, 10);
    m_peri_sigma ~ std_normal();
    zeta ~ normal(0, 20);

    // Likelihoods
    widths ~ normal(width_mu, width_sigma);
    alpha ~ normal(alpha_mu, alpha_sigma);
    rho_biomass ~ normal(rho_biomass_mu, rho_biomass_sigma);
    m_peri ~ normal(m_peri_mu, m_peri_sigma);

}


generated quantities {
    // PPcs
    vector<lower=0>[N_ppc] w_rep;
    vector<lower=0>[N_ppc] ell_rep;
    vector<lower=0>[N_ppc] vol_rep;
    vector<lower=0>[N_ppc] rho_peri_rep;

    for (i in 1:N_ppc) {
        w_rep[i] =  kappa * zeta .* (1 + phiRb_range[i] ./ beta_rp) / rho_biomass_mu;
        ell_rep[i] = alpha_mu * w_rep[i];
        vol_rep[i] = pi() * w_rep[i]^3 * (3 * alpha_mu - 1) / 12;
        rho_peri_rep[i] = m_peri_mu / (delta * pi() * alpha_mu * w_rep[i]^2);
    }

}