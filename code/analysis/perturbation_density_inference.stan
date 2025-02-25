data {
    // Dimensional parameters
    int<lower=1> N_prot; // Number of protein per cell measurements
    int<lower=1> N_size; // Number of cell size mesaurements
    int<lower=1> N_obs; // Number of experimental measurements on which to calculate empirics
    int<lower=1> N_ppc; // Number of posterior predictive checks for fits. 

    // Input data
    vector<lower=0>[N_prot] protein_per_cell;
    vector<lower=0>[N_prot] protein_per_cell_lam;
    vector<lower=0>[N_size] volume;
    vector<lower=0>[N_size] volume_lam;

    // Observed data
    vector<lower=0, upper=1>[N_obs] phi_rib;
    vector<lower=0, upper=1>[N_obs] phi_cyto;
    vector<lower=0, upper=1>[N_obs] phi_peri;
    vector<lower=0, upper=1>[N_obs] phi_mem;
    vector<lower=0>[N_obs] obs_surface_area;
    vector<lower=0>[N_obs] obs_volume;

    // PPC growth rate range
    vector<lower=0>[N_ppc] ppc_lam_range;

    // Constants
    real<lower=0> BETA_RIB; // Conversion factor from phi_rib to RNA-to-protein
    real<lower=0> W_PERI; // Periplasmic width in  microns
}

transformed data {
    vector[N_size] log_volume = log(volume);
}

parameters {
    // Regression parameters for exponential fits
    real<lower=0> volume_beta_0;
    real volume_beta_1; 

    // Homoscedastic errors in log transform
    real<lower=0> log_volume_sigma;
    real<lower=0> log_protein_sigma;

    // Inferred density
    real<lower=0> rho_prot_mu;

}

transformed parameters  {
    // Log transform the linear intercepts 
    real log_volume_beta_0 = log(volume_beta_0);

    // Compute the vector of densities given size realtionship and protein per cell 
    vector[N_prot] est_vol = log_volume_beta_0 + volume_beta_1 .* protein_per_cell_lam;
    
    // Compute the means for clarity
    vector[N_size] log_volume_mu = log_volume_beta_0 + volume_beta_1 .* volume_lam;
    vector[N_prot] prot_per_cell_mu = rho_prot_mu * exp(est_vol);
}

model { 
    // Set priors in linear space for regression. Parameter bounds chosen using
    //https://distribution-explorer.github.io/ with quantile setter mode. 
    // Percentile bounds are set as a comment in [2.5%, 97.5%]
    volume_beta_0 ~ inv_gamma(3.358, 0.7787); // [0.1 fL/cell, 1 fL/cell]
    volume_beta_1 ~ normal(0, 1);
    log_volume_sigma ~ normal(0, 1);
    rho_prot_mu ~ normal(300, 102); // [100 fg/fL, 500 fg/fL]
    log_protein_sigma ~ normal(0, 1);

    // Set likelihood
    log_volume ~ normal(log_volume_beta_0 + volume_beta_1 .* volume_lam, log_volume_sigma);
    log(protein_per_cell) ~ normal(log(prot_per_cell_mu), log_protein_sigma);
}

generated quantities { 
    // Posterior predictive checks from fit quantities
    vector[N_ppc] volume_ppc_mu;
    vector[N_ppc] prot_ppc_mu;
    vector[N_ppc] volume_ppc; 
    vector[N_ppc] protein_per_cell_ppc;

    // Empirical densities
    vector[N_obs] prot_per_cell;
    vector[N_obs] rna_per_cell;
    vector[N_obs] cyt_tot_per_cell;
    vector[N_obs] peri_prot_per_cell;
    vector[N_obs] mem_prot_per_cell;
    vector[N_obs] rho_cyt;
    vector[N_obs] rho_peri;
    vector[N_obs] sigma_mem;
    vector[N_obs] empirical_kappa;

    // Posterior predictive check
    for (i in 1:N_ppc) {
        volume_ppc_mu[i] = exp(log_volume_beta_0 + volume_beta_1 * ppc_lam_range[i]);
        prot_ppc_mu[i] = rho_prot_mu * volume_ppc_mu[i];
        volume_ppc[i] = exp(normal_rng(log(volume_ppc_mu[i]), log_volume_sigma));
        protein_per_cell_ppc[i] = exp(normal_rng(log(prot_ppc_mu[i]), log_protein_sigma));
    }

    // Generated quantities
    for (i in 1:N_obs) {
        prot_per_cell[i] = rho_prot_mu * obs_volume[i];
        rna_per_cell[i] = BETA_RIB * phi_rib[i] * prot_per_cell[i];
        cyt_tot_per_cell[i] = rna_per_cell[i] + phi_cyto[i] * prot_per_cell[i];
        peri_prot_per_cell[i] = phi_peri[i] * prot_per_cell[i];
        mem_prot_per_cell[i] = phi_mem[i] * prot_per_cell[i];
        rho_cyt[i] = cyt_tot_per_cell[i] / (obs_volume[i] - W_PERI * obs_surface_area[i]);
        rho_peri[i] = peri_prot_per_cell[i] / (W_PERI * obs_surface_area[i]);
        sigma_mem[i] = mem_prot_per_cell[i] / (2 * obs_surface_area[i]);
        empirical_kappa[i] = rho_cyt[i] / sigma_mem[i];
    }
}