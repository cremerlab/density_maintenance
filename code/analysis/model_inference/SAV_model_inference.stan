
data {
    // Test data
    int<lower=1> N_pred;
    vector<lower=0>[N_pred] pred_phiRb_range;
    vector<lower=0>[N_pred] pred_lam_range;

    //--------------------------------------------------------------------------
    // LITERATURE DATA
    //--------------------------------------------------------------------------

    // Literature size data.
    int<lower=1> N_size;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_size] surface_areas;
    vector<lower=0>[N_size] aspect_ratios;

    // Protein-per-cell data 
    int<lower=0> N_prot;
    vector<lower=0>[N_prot] prot_lam;
    vector<lower=0>[N_prot] prot_per_cell;

    // Drymass density data
    int<lower=0> N_drymass;
    vector<lower=0>[N_drymass] drymass_lam;
    vector<lower=0>[N_drymass] drymass_density;

    // Mass spec data
    int<lower=0> N_ms;
    vector<lower=0>[N_ms] ms_lam;
    vector<lower=0>[N_ms] phi_mem;
    vector<lower=0>[N_ms] phi_peri;

    // PhiRb data
    int<lower=1> N_phiRb;
    vector<lower=0>[N_phiRb] phiRb;
    vector<lower=0>[N_phiRb] phiRb_lam;
}

parameters {
    // Allocation parameters
    real<lower=0> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real<lower=0> phi_peri_sigma;
    real<lower=0> phiRb_intercept;
    real<lower=0> phiRb_slope;
    real<lower=0> phiRb_sigma;

    real<lower=0> m_peri;
    real<lower=0> drymass_mu;
    real<lower=0> drymass_sigma;

    // Shape parameters for calculating rho_mem from the mass spec data.
    real<lower=0> surface_area_slope;
    real<lower=0> surface_area_intercept;
    real<lower=0> surface_area_sigma;
    real<lower=0> alpha_zeroed_mu;
    real<lower=0> alpha_sigma;

    real<lower=0> rho_mem_mu;
    real log_prot_intercept;
    real<lower=0> log_prot_slope;
    real<lower=0> log_prot_sigma;

}

transformed parameters {
    vector<lower=0>[N_ms] ms_total_prot = exp(log_prot_intercept + log_prot_slope .* ms_lam);
    vector<lower=0>[N_ms] ms_phi_peri = m_peri ./ ms_total_prot;
    vector<lower=0>[N_ms] ms_surface_area = surface_area_intercept + surface_area_slope .* ms_lam; 
    real<lower=1> alpha_mu = alpha_zeroed_mu + 1;
    real<lower=0> kappa = drymass_mu / rho_mem_mu;

}

model { 

    // Membrane allocation and density
    rho_mem_mu ~ std_normal();
    phi_mem_mu ~ beta(2.5, 8.5);
    phi_mem_sigma ~ normal(0, 0.1);
    phi_mem ~ normal(phi_mem_mu, phi_mem_sigma);
    phi_mem ~ normal(rho_mem_mu .* 2 .* ms_surface_area ./ ms_total_prot, phi_mem_sigma);
     
    // Periplasmic allocation
    m_peri ~ normal(0, 10);
    phi_peri_sigma ~ normal(0, 0.1);
    phi_peri ~ normal(ms_phi_peri, phi_peri_sigma);

    // Ribosomal allocation
    phiRb_intercept ~ normal(0,  0.1);
    phiRb_slope ~ std_normal();
    phiRb ~ normal(phiRb_intercept + phiRb_slope .* phiRb_lam, phiRb_sigma);

    // Literature surface area inference
    surface_area_intercept ~ normal(0, 3);
    surface_area_slope ~ normal(0, 2);
    surface_area_sigma ~ std_normal();
    surface_areas ~ normal(surface_area_intercept + surface_area_slope .* size_lam, surface_area_sigma);

    // Literature aspect ratio inference
    alpha_zeroed_mu ~ normal(0, 3);
    alpha_sigma ~ std_normal();
    aspect_ratios ~ normal(alpha_mu, alpha_sigma);

    // Literature drymass inference
    drymass_mu ~ normal(300, 20);
    drymass_sigma ~ normal(0, 10);
    drymass_density ~ normal(drymass_mu, drymass_sigma);

    // Protein per cell
    log_prot_intercept ~ std_normal();
    log_prot_slope ~ normal(3, 2);
    log_prot_sigma ~ normal(0, 0.1);
    log(prot_per_cell) ~ normal(log_prot_intercept + log_prot_slope .* prot_lam, log_prot_sigma);

}

generated quantities {
    // Calculated theory parameters
    vector[N_pred] pred_phiRb = phiRb_intercept + phiRb_slope .* pred_lam_range;
    vector[N_pred] pred_lam = (pred_phiRb_range - phiRb_intercept)  ./ phiRb_slope;
    vector[N_pred] pred_phiRb_prot = exp(log_prot_intercept + log_prot_slope .* pred_lam);
    vector[N_pred] pred_lam_prot = exp(log_prot_intercept + log_prot_slope .* pred_lam_range);
    vector[N_pred] SAV_theory = (phi_mem_mu .* kappa) ./ (2 * (1 + (pred_phiRb_range / 0.4558) - phi_mem_mu - m_peri ./ pred_phiRb_prot)); 
    vector[N_pred] width_theory = (24 * alpha_mu / (3 * alpha_mu - 1)) .* (1 + (pred_phiRb_range/0.4558) - phi_mem_mu - m_peri ./pred_phiRb_prot) / (kappa * phi_mem_mu);
    vector[N_pred] width_pred = (24 * alpha_mu / (3 * alpha_mu - 1)) .* (1 + (pred_phiRb/0.4558) - phi_mem_mu - m_peri ./pred_lam_prot) / (kappa * phi_mem_mu);
    vector[N_pred] volume_pred = (pi() / 12) * width_pred^3 * (3 * alpha_mu  - 1);
    vector[N_pred] length_pred = alpha_mu * width_pred;
    vector[N_pred] aspect_ratio_pred;
    vector[N_pred] m_peri_pred;
    vector[N_pred] rho_mem_pred;
    vector[N_pred] phi_mem_pred;
    vector[N_pred] phi_peri_pred;
    for (i in 1:N_pred) {
        aspect_ratio_pred[i] = alpha_mu;
        m_peri_pred[i] = m_peri;
        rho_mem_pred[i] = rho_mem_mu;
        phi_mem_pred[i] = phi_mem_mu;
        phi_peri_pred[i] = m_peri ./ pred_phiRb_prot[i];
    } 

    // Empirical calculations from data
    vector[N_size] size_phiRb = phiRb_intercept + phiRb_slope .* size_lam;
    vector[N_ms] ms_m_peri = ms_total_prot .* phi_peri;
    vector[N_ms] ms_rho_mem = phi_mem .* ms_total_prot ./ (2 .* ms_surface_area);

    // // Modeled data posterior predictive checks
    vector[N_phiRb]  phiRb_ppc;
    vector[N_size] surface_area_ppc;
    real aspect_ratios_ppc = normal_rng(alpha_mu, alpha_sigma);
    vector[N_prot] prot_per_cell_ppc;
    vector[N_ms] phi_mem_mu_ppc;
    vector[N_ms] phi_mem_density_ppc;
    vector[N_ms] phi_peri_ppc;
    vector[N_drymass] drymass_density_ppc;

    for (i in 1:N_phiRb) {
        phiRb_ppc[i] = normal_rng(phiRb_intercept + phiRb_slope .* phiRb_lam[i], phiRb_sigma);
    }

    for (i in 1:N_size) {
        surface_area_ppc[i] = normal_rng(surface_area_intercept + surface_area_slope .* size_lam[i], surface_area_sigma);
    }

    for (i in 1:N_prot) { 
        prot_per_cell_ppc[i] = exp(normal_rng(log_prot_intercept + log_prot_slope * prot_lam[i], log_prot_sigma));
    }

    for (i in 1:N_ms) {
        phi_mem_density_ppc[i] = normal_rng(rho_mem_mu .* 2 .* ms_surface_area[i] ./ ms_total_prot[i], phi_mem_sigma);
        phi_mem_mu_ppc[i] = normal_rng(phi_mem_mu, phi_mem_sigma);
        phi_peri_ppc[i] = normal_rng(m_peri ./ exp(log_prot_intercept + log_prot_slope .* ms_lam[i]), phi_peri_sigma);
    }

    for (i in 1:N_drymass) {
        drymass_density_ppc[i] = normal_rng(drymass_mu, drymass_sigma);
    }

}
