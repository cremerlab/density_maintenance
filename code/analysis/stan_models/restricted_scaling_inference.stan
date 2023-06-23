data {
    // Dimensions
    int<lower=1> N_size;
    int<lower=1> N_alpha;
    int<lower=1> N_ms;
    int<lower=1> N_biomass;
    int<lower=1> N_ppc;

    // Observed   
    vector<lower=0>[N_size] width;
    vector<lower=0>[N_size] phiRb;
    vector<lower=1>[N_alpha] alpha;
    vector<lower=0>[N_ms] phi_mem;
    vector<lower=0>[N_ms] m_peri;
    vector<lower=0>[N_biomass] rho_biomass;

    // Calculated
    real<lower=0> beta_rp;
    real<lower=0> sa_prefactor;
    real<lower=0> delta;
    vector<lower=0>[N_ppc] phiRb_range;
}

// transformed data {
//     vector[N_biomass] rho_biomass_centered = (rho_biomass - mean(rho_biomass)) ./ sd(rho_biomass);
// }

parameters { 
    real<lower=0> width_sigma;
    real<lower=1> alpha_mu;
    real<lower=0> alpha_sigma;
    real<lower=0, upper=1> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    real<lower=0> rho_biomass_mu;
    real<lower=0> rho_biomass_sigma;
    real<lower=0> k;
}

transformed parameters {
    // vector<lower=0>[N_biomass] rho_biomass_mu = rho_biomass_centered .* sd(rho_biomass) + mean(rho_biomass); 
    vector<lower=0>[N_size] width_mu = (12 * sa_prefactor * alpha_mu) .* (1 + phiRb ./ beta_rp) ./(k * (3 * alpha_mu - 1) * phi_mem_mu);
}

model {
    // Priors
    m_peri_mu ~ std_normal();
    m_peri_sigma ~ std_normal();
    alpha_mu ~ std_normal();
    alpha_sigma ~ std_normal();
    width_sigma ~ std_normal();
    phi_mem_mu ~ normal(0, 0.1);
    phi_mem_sigma ~ normal(0, 0.1);
    rho_biomass_mu ~ normal(300, 100);
    rho_biomass_sigma ~ normal(0, 10);
    k ~ lognormal(6, 2);

    // Likelihoods
    m_peri ~ normal(m_peri_mu, m_peri_sigma);
    alpha ~ normal(alpha_mu, alpha_sigma);
    phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
    rho_biomass ~ normal(rho_biomass_mu, rho_biomass_sigma);
    width ~ normal(width_mu, width_sigma);
}

generated quantities {
    // Constant PPCs
    real<lower=1> alpha_rep = normal_rng(alpha_mu, alpha_sigma);
    real<lower=0> phi_mem_rep = normal_rng(phi_mem_mu, phi_mem_sigma);
    real<lower=0> rho_biomass_rep = normal_rng(rho_biomass_mu, rho_biomass_sigma);
    real<lower=0> m_peri_rep = normal_rng(m_peri_mu, m_peri_sigma);
    real<lower=0> rho_mem = rho_biomass_rep / k;

    // Simulated PPCs
    vector<lower=0>[N_ppc] width_rep;
    vector<lower=0>[N_ppc] ell_rep;
    vector<lower=0>[N_ppc] vol_rep;
    vector<lower=0>[N_ppc] rho_peri_rep;
    for (i in 1:N_ppc) {
        width_rep[i] = normal_rng((12 * sa_prefactor * alpha_mu) * (1 + phiRb_range[i] ./ beta_rp) ./ (k * (3 * alpha_mu - 1) * phi_mem_mu), width_sigma);
        ell_rep[i] = alpha_mu * width_rep[i];
        vol_rep[i] = pi() * width_rep[i]^3 * (3 * alpha_mu - 1) / 12;
        rho_peri_rep[i] = m_peri_mu / (pi() * alpha_mu * width_rep[i]^2 * delta);
    }
}