// Fits linear (and exponential) models to proteome scaling data 
data {

    // Dimensions
    int<lower=1> N_obs; // Number of measurements
    int<lower=1> N_pred; // Number of prediction points

    // Observed data
    vector<lower=0>[N_obs] growth_rate_hr; // Growth rate measurements
    vector<lower=0>[N_obs] phi_rib; // Ribosomal fraction measurements
    vector<lower=0>[N_obs] phi_peri; // Periplasmic fraction measurements
    vector<lower=0>[N_obs] phi_mem; // Membrane fraction measurements

    // Predition data
    vector<lower=0>[N_pred] growth_rate_hr_pred; // Growth rate predictions
}

transformed data {
    // Log transform of the phi_peri data so we can fit a simple linear model
    vector[N] log_phi_peri;
}

parameters { 
    real beta_1_rib; // Slope of the ribosomal scaling relationship
    real beta_0_rib; // Intercept of the ribosomal scaling relationship
    real<lower=0> phi_rib_sigma; // Measurement error of the ribosomal scaling relationship

    real beta_1_peri; // Slope of the log-periplasmic scaling relationship
    real beta_0_peri; // Intercept of the periplasmic scaling relationship
    real<lower=0> phi_peri_sigma; // Measurement error of the log-periplasmic scaling relationship

    real beta_1_mem; // Slope of the membrane scaling relationship
    real beta_0_mem; // Intercept of the membrane scaling relationship
    real<lower=0> phi_mem_sigma; // Measurement error of the membrane scaling relationship
}

transformed parameters {
    real log_beta_0_peri = log(beta_0_peri); // For more intuitive setting of priors
}

model {
    // Define the priors using distribution-explorer.github.io 
    beta_1_rib ~ normal(0, 1);     
    beta_0_rib ~ beta(3.242, 71.08); // Quantiles [2.5%, 97.5%] = [0.01, 0.1]
    phi_rib_sigma ~ normal(0, 0.1); 
    beta_1_peri ~ normal(0, 1); 
    beta_0_peri ~ beta(2.037, 17.74); // Quantiles [2.5%, 97.5%] = [0.01, 0.25]
    phi_peri_sigma ~ normal(0, 1);
    beta_1_mem ~ normal(0, 1);
    beta_0_mem ~ beta(2.037, 17.74); // Quantiles [2.5%, 97.5%] = [0.01, 0.25]
    phi_mem_sigma ~ normal(0, 0.1);

    // Define the likelihoods
    phi_rib ~ normal(beta_0_rib + beta_1_rib * growth_rate_hr, phi_rib_sigma);
    log_phi_peri ~ normal(log_beta_0_peri + beta_1_peri * growth_rate_hr, phi_peri_sigma);
    phi_mem ~ normal(beta_0_mem + beta_1_mem * growth_rate_hr, phi_mem_sigma);
}

generated quantities {
    // Generate posterior preditive checks for the scaling relationships
    vector[N_pred] phi_rib_ppc;
    vector[N_pred] phi_peri_ppc;
    vector[N_pred] phi_mem_ppc;

    for (i in 1:N_pred) {
        phi_rib_ppc[i] = normal_rng(beta_0_rib + beta_1_rib * growth_rate_hr_pred[i], phi_rib_sigma);
        phi_peri_ppc[i] = exp(normal_rng(log_beta_0_peri + beta_1_peri * growth_rate_hr_pred[i], phi_peri_sigma));
        phi_mem_ppc[i] = normal_rng(beta_0_mem + beta_1_mem * growth_rate_hr_pred[i], phi_mem_sigma);
    }
}
