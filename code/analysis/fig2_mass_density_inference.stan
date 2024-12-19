data { 
    // Input dimensions
    int<lower=1> N_prot; // Number of protein per cell measurements
    int<lower=1> N_size; // Number of size measurement
    int<lower=1> N_ms; // Number of mass spec measurements
    int<lower=1> N_obs; // Number of experimental measurements made for this work
    int<lower=1> N_fit; // Number of data points for fit generation

    // Input independent data for inference
    vector<lower=0>[N_prot] prot_per_cell;
    vector<lower=0>[N_size] surface_area;
    vector<lower=0>[N_size] volume;

    // Input dependent data for inference
    vector<lower=0>[N_prot] prot_per_cell_lam;
    vector<lower=0>[N_size] size_lam;

    // Input data for generative modeling
    vector<lower=0>[N_ms] phi_cyto;
    vector<lower=0>[N_ms] phi_mem;
    vector<lower=0>[N_ms] phi_peri;
    vector<lower=0>[N_ms] phi_rib;
    vector<lower=0>[N_ms] ms_lam;

    // Input experimental data for density calculations
    vector<lower=0>[N_obs] obs_volume;
    vector<lower=0>[N_obs] obs_surface_area;
    vector<lower=0>[N_obs] obs_phi_cyto;
    vector<lower=0>[N_obs] obs_phi_mem;
    vector<lower=0>[N_obs] obs_phi_peri;
    vector<lower=0>[N_obs] obs_phi_rib;
    vector<lower=0>[N_obs] obs_lam;

    // Input dims for continuous generative fits
    vector<lower=0>[N_fit] fit_lam;

    // Constants for calculation
    real<lower=0> BETA; // For conversion from mass frac to RNA/Protein
    real<lower=0> rRNA_FRAC; // Fraction of all RNA that is ribosomal.
    real<lower=0> W_PERI; // Width of the periplasm for calculation of peri volume.
}

transformed data {
    // Log transformations for easier inference
    vector[N_prot] log_prot_per_cell = log(prot_per_cell);
    vector[N_size] log_volume = log(volume);
    vector[N_size] log_surface_area = log(surface_area);
}

parameters { 
    // Regression parameters for exponential fits.
    real<lower=0> prot_per_cell_beta_0;  
    real prot_per_cell_beta_1;
    real<lower=0> volume_beta_0;
    real volume_beta_1;
    real<lower=0> surface_area_beta_0;
    real surface_area_beta_1;

    // Homoscedastic errors in the log transform
    real<lower=0> log_prot_per_cell_sigma;
    real<lower=0> log_volume_sigma;
    real<lower=0> log_surface_area_sigma;

}

transformed parameters {
    // Log transform the linear intercepts
    real log_prot_per_cell_beta_0 = log(prot_per_cell_beta_0);
    real log_volume_beta_0 = log(volume_beta_0);
    real log_surface_area_beta_0 = log(surface_area_beta_0);

    // Compute the means for clarity
    vector[N_prot] log_prot_per_cell_mu = log_prot_per_cell_beta_0 + prot_per_cell_beta_1 .* prot_per_cell_lam;
    vector[N_size] log_volume_mu = log_volume_beta_0 + volume_beta_1 .* size_lam;
    vector[N_size] log_surface_area_mu = log_surface_area_beta_0 + surface_area_beta_1 .* size_lam;
}

model {
    // Set priors in linear space for regression. Parameter bounds chosen using
    //https://distribution-explorer.github.io/ with quantile setter mode. 
    // Percentile bounds are set as a comment in [2.5%, 97.5%]
    prot_per_cell_beta_0 ~ inv_gamma(13.21, 1062); // [50 fg/cell , 150 fg/cell]
    prot_per_cell_beta_1 ~ normal(0, 1);
    volume_beta_0 ~ inv_gamma(2.153, 0.5836); // [0.1 fL, 2 fL]
    volume_beta_1 ~ normal(0, 1);
    surface_area_beta_0 ~ inv_gamma(6.407, 12.24); // [1 µm^2, 5 µm^2]
    surface_area_beta_1 ~ normal(0, 1);
    log_prot_per_cell ~ normal(0, 1);
    log_volume_sigma ~ normal(0, 1);
    log_surface_area_sigma ~ normal(0, 1);

    // Set likelihoods for inference
    log_prot_per_cell ~ normal(log_prot_per_cell_mu, log_prot_per_cell_sigma);
    log_volume ~ normal(log_volume_mu, log_volume_sigma);
    log_surface_area ~ normal(log_surface_area_mu, log_surface_area_sigma);
}

generated quantities {
    // Define empirical quantities to compute.
    vector[N_ms] emp_M_cyto;
    vector[N_ms] emp_M_peri;
    vector[N_ms] emp_M_mem;
    vector[N_ms] emp_rho_cyto;
    vector[N_ms] emp_rho_peri;
    vector[N_ms] emp_sigma_mem;
    vector[N_ms] emp_prot_per_cell_ms;
    vector[N_ms] emp_surface_area_ms;
    vector[N_ms] emp_volume_ms;
    vector[N_ms] emp_rho_rib;
    vector[N_ms] emp_rho_rrna;
    vector[N_ms] emp_rho_rna;
    vector[N_ms] emp_rho_biomass;


    // Define observed quantities for our mass spec data
    vector[N_obs] obs_M_cyto;
    vector[N_obs] obs_M_peri;
    vector[N_obs] obs_M_mem;
    vector[N_obs] obs_prot_per_cell;
    vector[N_obs] obs_rho_cyto;
    vector[N_obs] obs_rho_peri;
    vector[N_obs] obs_sigma_mem;
    vector[N_obs] obs_rho_rib;
    vector[N_obs] obs_rho_rrna;
    vector[N_obs] obs_rho_rna;
    vector[N_obs] obs_rho_biomass;

    // Define posterior predictive checks on fit range.
    vector[N_fit] prot_per_cell_ppc;
    vector[N_fit] volume_ppc;
    vector[N_fit] surface_area_ppc;
    
    // Compute empirical quantities
    for (i in 1:N_ms) {
        emp_prot_per_cell_ms[i] = prot_per_cell_beta_0 * exp(prot_per_cell_beta_1 * ms_lam[i]);   
        emp_surface_area_ms[i] = surface_area_beta_0 * exp(surface_area_beta_1 * ms_lam[i]); 
        emp_volume_ms[i] = volume_beta_0 * exp(volume_beta_1 * ms_lam[i]);
        emp_M_cyto[i] = phi_cyto[i] * emp_prot_per_cell_ms[i];
        emp_M_peri[i] = phi_peri[i] * emp_prot_per_cell_ms[i];
        emp_M_mem[i] = phi_mem[i] * emp_prot_per_cell_ms[i];
        emp_rho_cyto[i] = emp_M_cyto[i] / (emp_volume_ms[i] - emp_surface_area_ms[i] * W_PERI);
        emp_rho_peri[i] = emp_M_peri[i] / (emp_surface_area_ms[i] * W_PERI);
        emp_sigma_mem[i] = emp_M_mem[i] / (2 * emp_surface_area_ms[i]);
        emp_rho_rib[i] = phi_rib[i] * emp_prot_per_cell_ms[i] / (emp_volume_ms[i] - emp_surface_area_ms[i] * W_PERI);
        emp_rho_rrna[i] = BETA * emp_rho_rib[i];
        emp_rho_rna[i] = emp_rho_rrna[i] / rRNA_FRAC;
        emp_rho_biomass[i] =  emp_rho_rna[i] + emp_rho_cyto[i];
        }

    // Compute observed quantities 
    for (i in 1:N_obs) {
        obs_prot_per_cell[i] = prot_per_cell_beta_0 * exp(prot_per_cell_beta_1 * obs_lam[i]);
        obs_M_cyto[i] = obs_phi_cyto[i] * obs_prot_per_cell[i];
        obs_M_peri[i] = obs_phi_peri[i] * obs_prot_per_cell[i];
        obs_M_mem[i] = obs_phi_mem[i] * obs_prot_per_cell[i];
        obs_rho_cyto[i] = obs_M_cyto[i] / (obs_volume[i] - obs_surface_area[i] * W_PERI);
        obs_rho_peri[i] = obs_M_peri[i] / (obs_surface_area[i] * W_PERI);
        obs_sigma_mem[i] = obs_M_mem[i] / (2 * obs_surface_area[i]);
        obs_rho_rib[i] = obs_phi_rib[i] * obs_prot_per_cell[i] / (obs_volume[i] - obs_surface_area[i] * W_PERI);
        obs_rho_rrna[i] = BETA * obs_rho_rib[i];
        obs_rho_rna[i] = obs_rho_rrna[i] / rRNA_FRAC;
        obs_rho_biomass[i] = obs_rho_rna[i] + obs_rho_cyto[i];
    }

    // Compute posterior predictive checks
    for (i in 1:N_fit) {
        prot_per_cell_ppc[i] = exp(normal_rng(log_prot_per_cell_beta_0 + prot_per_cell_beta_1 .* fit_lam[i], log_prot_per_cell_sigma));
        surface_area_ppc[i] = exp(normal_rng(log_surface_area_beta_0 + surface_area_beta_1 .* fit_lam[i], log_surface_area_sigma));
        volume_ppc[i] = exp(normal_rng(log_volume_beta_0 + volume_beta_1 * fit_lam[i], log_volume_sigma));
    }

}

