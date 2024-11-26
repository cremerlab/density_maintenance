data { 
    // Input dimensions
    int<lower=1> N_prot; // Number of protein per cell measurements
    int<lower=1> N_size; // Number of size measurement
    int<lower=1> N_ms; // Number of mass spec measurements
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
    vector<lower=0>[N_ms] ms_lam;

    // Input dims for continuous generative fits
    vector<lower=0>[N_fit] fit_lam;

    // Constants for calculation
    real<lower=0> W_PERI;
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
    vector[N_ms] M_cyto;
    vector[N_ms] M_peri;
    vector[N_ms] M_mem;
    vector[N_ms] rho_cyto;
    vector[N_ms] rho_peri;
    vector[N_ms] sigma_mem;
    vector[N_ms] prot_per_cell_ms;
    vector[N_ms] surface_area_ms;
    vector[N_ms] volume_ms;

    // Define posterior predictive checks on fit range.
    vector[N_fit] prot_per_cell_ppc;
    vector[N_fit] volume_ppc;
    vector[N_fit] surface_area_ppc;
    
    // Compute empirical quantities
    for (i in 1:N_ms) {
        prot_per_cell_ms[i] = prot_per_cell_beta_0 * exp(prot_per_cell_beta_1 * ms_lam[i]);   
        surface_area_ms[i] = surface_area_beta_0 * exp(surface_area_beta_1 * ms_lam[i]); 
        volume_ms[i] = volume_beta_0 * exp(volume_beta_1 * ms_lam[i]);
        M_cyto[i] = phi_cyto[i] * prot_per_cell_ms[i];
        M_peri[i] = phi_peri[i] * prot_per_cell_ms[i];
        M_mem[i] = phi_mem[i] * prot_per_cell_ms[i];
        rho_cyto[i] = M_cyto[i] / (volume[i] - surface_area[i] * W_PERI);
        rho_peri[i] = M_peri[i] / (surface_area[i] * W_PERI);
        sigma_mem[i] = M_mem[i] / (2 * surface_area[i]);
        }

    // Compute posterior predictive checks
    for (i in 1:N_fit) {
        prot_per_cell_ppc[i] = exp(normal_rng(log_prot_per_cell_beta_0 + prot_per_cell_beta_1 .* fit_lam[i], log_prot_per_cell_sigma));
        surface_area_ppc[i] = exp(normal_rng(log_surface_area_beta_0 + surface_area_beta_1 .* fit_lam[i], log_surface_area_sigma));
        volume_ppc[i] = exp(normal_rng(log_volume_beta_0 + volume_beta_1 * fit_lam[i], log_volume_sigma));
    }
}

