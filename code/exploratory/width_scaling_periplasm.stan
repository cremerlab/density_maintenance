data {
    int<lower=0> N_prot;
    int<lower=0> N_phiRb;
    int<lower=0> N_ms;
    int<lower=0> N_biomass;
    int<lower=0> N_dna;
    int<lower=0> N_size;
    int<lower=0> N_ppc;

    vector<lower=0>[N_prot] prot_per_cell;
    vector<lower=0>[N_prot] prot_per_cell_lam;
    vector<lower=0>[N_phiRb] phiRb;
    vector<lower=0>[N_phiRb] phiRb_lam;
    vector<lower=0>[N_dna] dna_to_protein;
    vector<lower=0>[N_dna] dna_lam;
    vector<lower=0>[N_ms] phi_mem;
    vector<lower=0>[N_ms] rho_mem;
    vector<lower=0>[N_ms] m_peri;
    vector<lower=0>[N_biomass] rho_biomass;
    vector<lower=0>[N_size] aspect_ratio;
    vector<lower=0>[N_size] width;
    vector<lower=0>[N_size] vol;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_ppc] lam_range;
    vector<lower=0>[N_ppc] phiRb_range;

}

transformed data {
    vector[N_biomass] rho_biomass_centered = (rho_biomass - mean(rho_biomass)) / sd(rho_biomass);
}

parameters {
    // Constants
    real<lower=0, upper=1>phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=1> alpha_mu;
    real<lower=0> alpha_sigma;
    real<lower=0> rho_mem_mu;
    real<lower=0> rho_mem_sigma;
    real rho_biomass_centered_mu;
    real<lower=0> rho_biomass_centered_sigma;
    real<lower=0> m_peri_mu;
    real<lower=0> m_peri_sigma;
    real<lower=0> dna_mu;
    real<lower=0> dna_sigma;

    // Relationships
    real<lower=0> prot_per_cell_slope;
    real prot_per_cell_intercept;
    real<lower=0> prot_per_cell_sigma;

    real<lower=0> width_slope;
    real<lower=0> width_intercept;
    real<lower=0> width_sigma;

    real<lower=0> phiRb_slope;
    real<lower=0> phiRb_intercept;
    real<lower=0> phiRb_sigma;

    real<lower=0> vol_slope;
    real vol_intercept;
    real<lower=0> vol_sigma;


}

transformed parameters {
    real rho_biomass_mu = rho_biomass_centered_mu * sd(rho_biomass) + mean(rho_biomass);
    real kappa = rho_biomass_mu / rho_mem_mu;

}

model { 

    // Membrane density
    phi_mem_mu ~ std_normal();
    phi_mem_sigma ~ std_normal();
    phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);

    // Biomass density
    rho_biomass_centered_mu ~ std_normal();
    rho_biomass_centered_sigma ~ std_normal();
    rho_biomass_centered ~ normal(rho_biomass_centered_mu, rho_biomass_centered_sigma);

    // Membrane density
    rho_mem_mu ~ std_normal();
    rho_mem_sigma ~ std_normal();
    rho_mem ~ normal(rho_mem_mu, rho_mem_sigma);

    // Aspect ratio
    alpha_mu ~ normal(1, 1);
    alpha_sigma ~ std_normal();
    aspect_ratio ~ normal(alpha_mu, alpha_sigma);

    // Periplasmic mass
    m_peri_mu ~ normal(0, 10);
    m_peri_sigma ~ std_normal();
    m_peri ~ normal(m_peri_mu, m_peri_sigma);

    // DNA to protein ratio
    dna_mu ~ std_normal();
    dna_sigma ~ std_normal();
    dna_to_protein ~ normal(dna_mu, dna_sigma);

    // Protein per cell
    prot_per_cell_slope ~ std_normal();
    prot_per_cell_intercept ~ std_normal();
    prot_per_cell_sigma ~ std_normal();
    log(prot_per_cell) ~ normal(prot_per_cell_intercept + prot_per_cell_slope .* prot_per_cell_lam, prot_per_cell_sigma);

    // Width scaling
    width_intercept ~ std_normal();
    width_slope ~ std_normal();
    width_sigma ~ std_normal();
    width ~ normal(width_intercept + width_slope .* size_lam, width_sigma);

    // Volume scaling
    vol_intercept ~ std_normal();
    vol_slope ~ std_normal();
    vol_sigma ~ std_normal();
    log(vol) ~ normal(vol_intercept + vol_slope .* size_lam, vol_sigma);

    // R-line
    phiRb_intercept ~ std_normal();
    phiRb_slope ~ std_normal();
    phiRb_sigma ~ std_normal();
    phiRb ~ normal(phiRb_intercept + phiRb_slope * phiRb_lam, phiRb_sigma);
}

generated quantities {
    vector[N_ppc] phiRb_rep;
    vector[N_ppc] phiRb_sim;
    vector[N_ppc] width_rep;
    vector[N_ppc] width_sim;
    vector[N_ppc] prot_per_cell_rep;
    vector[N_ppc] prot_per_cell_sim;
    vector[N_ppc] vol_pred_sim;
    vector[N_ppc] vol_pred_rep;
    vector[N_ppc] width_pred_sim;
    vector[N_ppc] width_pred_rep;
    vector[N_ppc] rho_peri_rep;
    vector[N_ppc] rho_peri_sim;
    vector[N_ppc] phi_peri_sim;
    vector[N_ppc] phi_peri_rep;
    real alpha_rep = normal_rng(alpha_mu, alpha_sigma);
    real m_peri_rep = normal_rng(m_peri_mu, m_peri_sigma);
    real rho_mem_rep = normal_rng(rho_mem_mu, rho_mem_sigma);
    real rho_biomass_rep = mean(rho_biomass) + sd(rho_biomass) * normal_rng(rho_biomass_centered_mu, rho_biomass_centered_sigma);
    real phi_mem_rep = normal_rng(phi_mem_mu, phi_mem_sigma);
    real dna_to_protein_rep = normal_rng(dna_mu, dna_sigma);

    for (i in 1:N_ppc) {
        phiRb_sim[i] = phiRb_intercept + phiRb_slope * lam_range[i];
        phiRb_rep[i] = normal_rng(phiRb_sim[i], phiRb_sigma);
        width_sim[i] = width_intercept + width_slope * lam_range[i];
        width_rep[i] = normal_rng(width_sim[i], width_sigma);
        prot_per_cell_sim[i] = exp(prot_per_cell_intercept + prot_per_cell_slope * lam_range[i]);
        prot_per_cell_rep[i] = exp(normal_rng(log(prot_per_cell_sim[i]), prot_per_cell_sigma));
        phi_peri_sim[i] = m_peri_mu / prot_per_cell_sim[i]; 
        phi_peri_rep[i] = normal_rng(phi_peri_sim[i], phi_mem_sigma);
        width_pred_sim[i] = 0.0245 * alpha_mu + ((24 * alpha_mu / (3 * alpha_mu - 1)) * (1/(kappa * phi_mem_mu)) * (1 + (phiRb_range[i]/0.4558) - phi_mem_mu - phi_peri_sim[i] - dna_mu));
        width_pred_rep[i] = normal_rng(width_pred_sim[i], width_sigma);
        rho_peri_sim[i] = m_peri_mu / (pi() * alpha_mu * width_pred_sim[i]^2 * 0.0245);
        rho_peri_rep[i] = normal_rng(rho_peri_sim[i], rho_mem_sigma);
        vol_pred_sim[i] = (pi()/12) * width_pred_sim[i]^3 * (3 * alpha_mu - 1);
        vol_pred_rep[i] = exp(normal_rng(log(vol_pred_sim[i]), vol_sigma));
    }
    
}