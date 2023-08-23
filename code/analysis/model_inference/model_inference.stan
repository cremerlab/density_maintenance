//-----------------------------------------------------------------------------
// Model Inference
// ==================
// This model infers the central model parameter constants given literature
// measurements. This version propagates uncertainity through empirical 
// transformation of mass spectrometry data to obtain membrane density and 
// periplasmic protein mass. Uncertainty is also propagated for inference of 
// average cell width from the growth rate from measurements of the ribosomal 
// allocation.
//-----------------------------------------------------------------------------
data { 
    // Literature size data.
    int<lower=1> N_size;
    vector<lower=0>[N_size] size_lam;
    vector<lower=0>[N_size] widths;
    vector<lower=0>[N_size] aspect_ratio;

    // Protein-per-cell data 
    int<lower=0> N_prot;
    vector<lower=0>[N_prot] prot_lam;
    vector<lower=0>[N_prot] prot_per_cell;

    // Drymass density data
    int<lower=0> N_drymass;
    vector<lower=0>[N_drymass] drymass_lam;
    vector<lower=0>[N_prot] drymass_density;

    // Mass spec data
    int<lower=0> N_ms;
    vector<lower=0>[N_ms] ms_lam;
    vector<lower=0>[N_ms] phi_mem;
    vector<lower=0>[N_ms] phi_peri;

    // PhiRb data
    int<lower=1> N_phi;
    vector<lower=0>[N_phi] phi_lam;

    // Prediction data
    int<lower=1> N_pred;
    vector<lower=0> pred_lam;

}

parameters {
    // Inference of growth-rate dependent cell size inference
    real<lower=0> width_slope; 
    real width_intercept; 
    real<lower=0> width_sigma; 
    real<lower=0> alpha_zeroed_mu;
    real<lower=0> alpha_zeroed_sigma;

    // Protein per cell inference
    real log_prot_slope;
    real log_prot_intercept; 
    real log_prot_sigma;

    // Constant parameters of model.
    real<lower=1, upper=1> phi_mem_mu;
    real<lower=0> phi_mem_sigma;
    real log_rho_bio;
    real<lower=0> rho_bio_sigma;
    real<lower=0> m_peri;
    real<lower=0> phi_peri_sigma;
    real log_kappa; 
}

transformed parameters {
    // Empirical parameters
    vector<lower=0>[N_ms] ms_prot = exp(log_prot_intercept + log_prot_slope .* ms_lam);
    vector<lower=0>[N_ms] ms_width = width_intercept + width_slope .* ms_lam;
    vector<lower=0>[N_ms] ms_length = length_intercept + length_slope .* ms_lam;
    vector<lower=0>[N_ms] ms_surface_area = pi() .* ms_length .* ms_width;
    vector<lower=0>[N_ms] m_mem = phi_mem .* ms_prot; 
    vector<lower=0>[N_ms] m_peri = phi_peri .* ms_prot;
    vector<lower=0>[N_ms] rho_mem = m_mem ./ (2 .* ms_surface_area);
    vector<lower=0>[N_ms] rho_peri = m_peri ./ (ms_surface_area .* 0.0246);
    vector<lower=0>[N_phi] phi_width = width_intercept + width_slope .* phi_lam;

    // Undoing log transforms
    real<lower=0> rho_bio = exp(log_rho_bio);
    real<lower=0> kappa = exp(log_kappa);
    real<lower=1> alpha_mu= 1 + alpha_zeroed;

}

model  {
    // Size information
    width_slope ~ std_normal();
    width_intercept ~ std_normal();
    width_sigma ~ std_normal();
    alpha_zeroed_mu ~ normal(0, 2);
    alpha_zeroed_sigma ~ std_normal();
    widths ~ normal(width_intercept + width_slope .* size_lam, width_sigma);
    aspect_ratio - 1 ~ normal(alpha_zeroed_mu, alpha_zeroed_sigma);
    
    // Protein per cell
    log_prot_intercept ~ std_normal();
    log_prot_slope ~ std_normal();
    log_prot_sigma ~ std_normal(); 
    log(prot_per_cell) ~ normal(log_prot_intercept + log_prot_slope .* prot_lam, log_prot_sigma);

    // Drymass density
    log_rho_bio ~ normal(0, 5);
    rho_bio_sigma ~ normal(0, 10);
    drymass_density ~ normal_rng(rho_bio, rho_bio_sigma);

    // Membrane protein allocation inference
    log_kappa ~ normal(0, 3);
    phi_mem_mu ~ beta(2.25, 8);
    phi_mem_sigma ~ normal(0, 0.1);
    phi_mem ~ normal_rng(phi_mem_mu, phi_mem_sigma);
    phi_mem ~ normal_rng(rho_bio .* 2 .* pi() * alpha_mu .* (width_intercept + width_slope .* ms_lam) ./ (kappa * exp(log_prot_intercept + log_prot_slope .* ms_lam)), phi_mem_sigma);

    // Periplasm protein allocation inference  
    m_peri ~ normal(0, 5); 
    phi_peri ~ normal_rng(m_peri ./ exp(log_prot_intercept + log_prot_slope .* ms_lam), phi_peri_sigma);
}

generated quantities { 
    // Posterior predictive checks
    vector<lower=0>[N_size] width_ppc;
    vector<lower=0>[N_size] alpha_ppc;
    vector<lower=0>[N_prot] prot_per_cell_ppc;

    // Model predictions
    vector<lower=0>[N_pred] rho_mem_pred;
    vector<lower=0>[N_pred] width_pred;

    // Empirical data transformations
    vector<lower=0>[N_ms] rho_mem_ms;
    vector<lower=0>[N_ms] m_peri_ms;
    vector<lower=0>[N_phi] width_phiRb;

    for (i in 1:N_size) {
        width_ppc[i] = normal_rng(width_intercept + width_slope * size_lam[i], width_sigma);
        alpha_ppc[i] = 1 + normal_rng(alpha_zeroed_mu, alpha_zeroed_sigma);
    }

    for (i in 1:N_prot) {
        prot_per_cell_ppc[i] = exp(normal_rng(log_prot_intercept + log_prot_slope * prot_lam[i], log_prot_sigma));
    }

    for (i in 1:N_ms) {
        
    }

}