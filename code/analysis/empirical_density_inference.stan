data { 
     // Input dimensions
    int<lower=1> N_rp; // Number of protein per cell measurements from our experiments.
    int<lower=1> N_obs; // Number of mass spec measurements from our experiments.
    int<lower=1> N_fit; // Number of data points for fit generation   

    // Input the independent data for inference
    vector<lower=0>[N_rp] prot_per_cell; // Protein per cell measurements
    vector<lower=0>[N_rp] rna_per_cell; // RNA per cell measurements 
    vector<lower=0>[N_obs] surface_area; // Surface area measurements
    vector<lower=0>[N_obs] volume; // Volume measurements

    // Define the dependent data
    vector<lower=0>[N_rp] rp_lam; // Protein and RNA per cell growth rate
    vector<lower=0>[N_obs] lam; //  Growth rate
    
    // Define input data for generative modeling
    vector<lower=0>[N_obs] phi_cyto; // Cytoplasmic mass fraction measurements
    vector<lower=0>[N_obs] phi_mem; // Membrane mass fraction measurements
    vector<lower=0>[N_obs] phi_peri; // Periplasmic mass fraction measurements

    // Input data for generative fit generation
    vector<lower=0>[N_fit] fit_lam; // Growth rate measurements for fit generation
 
    // Constants
    real<lower=0> W_PERI; // Periplasmic width
}

transformed data {
    // Log transformations for easier inference
    vector[N_rp] log_prot_per_cell = log(prot_per_cell);
    vector[N_rp] log_rna_per_cell = log(rna_per_cell);

}

parameters {
    // Regression parameters for exponential fits.
    real<lower=0> prot_per_cell_beta_0;  
    real prot_per_cell_beta_1;
    real<lower=0> rna_per_cell_beta_0;
    real rna_per_cell_beta_1;

    // Homoscedastic errors in the log transform
    real<lower=0> log_prot_per_cell_sigma;
    real<lower=0> log_rna_per_cell_sigma;

}

transformed parameters { 
    // Log transform the linear intercepts
    real log_prot_per_cell_beta_0 = log(prot_per_cell_beta_0);
    real log_rna_per_cell_beta_0 = log(rna_per_cell_beta_0);

    // Compute the means for clarity
    vector[N_rp] log_prot_per_cell_mu = log_prot_per_cell_beta_0 + prot_per_cell_beta_1 .* rp_lam;
    vector[N_rp] log_rna_per_cell_mu = log_rna_per_cell_beta_0 + rna_per_cell_beta_1 .* rp_lam;

}

model { 
    // Set priors in linear space for regression. Parameter bounds chosen using
    //https://distribution-explorer.github.io/ with quantile setter mode. 
    // Percentile bounds are set as a comment in [2.5%, 97.5%]
    prot_per_cell_beta_0 ~ inv_gamma(13.21, 1062); // [50 fg/cell , 150 fg/cell]
    prot_per_cell_beta_1 ~ normal(0, 1);
    rna_per_cell_beta_0 ~ inv_gamma(3.358, 77.87);  // [10 fg / cell, 100 fg / cell]
    rna_per_cell_beta_1 ~ normal(0, 1);
    log_prot_per_cell_sigma ~ normal(0, 1);
    log_rna_per_cell_sigma ~ normal(0, 1);


    // Set likelihoods for inference
    log_prot_per_cell ~ normal(log_prot_per_cell_mu, log_prot_per_cell_sigma);
    log_rna_per_cell ~ normal(log_rna_per_cell_mu, log_rna_per_cell_sigma);
}

generated quantities {
    // Posterior predictive checks for fit quantities
    vector[N_fit] prot_per_cell_ppc;
    vector[N_fit] rna_per_cell_ppc;

    // Empirical quantities from our measurements
    vector[N_obs] tot_prot_per_cell;
    vector[N_obs] cyt_prot_per_cell;
    vector[N_obs] cyt_rna_per_cell;
    vector[N_obs] cyt_tot_per_cell;
    vector[N_obs] peri_prot_per_cell;
    vector[N_obs] mem_prot_per_cell;
    vector[N_obs] rho_cyt_prot;
    vector[N_obs] rho_cyt_rna;
    vector[N_obs] rho_cyt_tot;
    vector[N_obs] rho_peri;
    vector[N_obs] sigma_mem;
    vector[N_obs] empirical_kappa;

    // Compute the posterior predictive checks.
    for (i in 1:N_fit) { 
        prot_per_cell_ppc[i] = exp(normal_rng(log_prot_per_cell_beta_0 + prot_per_cell_beta_1 * fit_lam[i], log_prot_per_cell_sigma)); 
        rna_per_cell_ppc[i] = exp(normal_rng(log_rna_per_cell_beta_0 + rna_per_cell_beta_1 * fit_lam[i], log_rna_per_cell_sigma));
    }

    // Compute the empirical quantities. 
    for (i in 1:N_obs) {
        tot_prot_per_cell[i] = exp(log_prot_per_cell_beta_0 + prot_per_cell_beta_1 * lam[i]);
        cyt_rna_per_cell[i] = exp(log_rna_per_cell_beta_0  + rna_per_cell_beta_1 * lam[i]);
        cyt_prot_per_cell[i] = phi_cyto[i] * tot_prot_per_cell[i];
        cyt_tot_per_cell[i] = cyt_prot_per_cell[i] + cyt_rna_per_cell[i]; 
        peri_prot_per_cell[i] = phi_peri[i] * tot_prot_per_cell[i];
        mem_prot_per_cell[i] = phi_mem[i] * tot_prot_per_cell[i];
        rho_cyt_prot[i] = cyt_prot_per_cell[i] / (volume[i] - W_PERI * surface_area[i]);
        rho_cyt_rna[i] = cyt_rna_per_cell[i] / (volume[i] - W_PERI * surface_area[i]);
        rho_cyt_tot[i] = (cyt_prot_per_cell[i] + cyt_rna_per_cell[i]) / (volume[i] - W_PERI * surface_area[i]);
        rho_peri[i] = peri_prot_per_cell[i] / (W_PERI * surface_area[i]);
        sigma_mem[i] = mem_prot_per_cell[i] / (2 * surface_area[i]);
        empirical_kappa[i] = rho_cyt_tot[i] / sigma_mem[i];
    }      
}