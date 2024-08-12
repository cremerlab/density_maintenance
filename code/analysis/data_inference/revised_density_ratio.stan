data {
    // Dimensional information for fitting 
    int<lower=1> N_prot; // protein per cell data
    int<lower=1> N_SA; // surface-to-volume measurements
    int<lower=1> N_vol; // volume measurements
    int<lower=1> N_DNA; // DNA-to-protein ratio measurements

    int<lower=1> N_scan;
    vector<lower=0>[N_scan] scan_lam;
    // Dimensional information for calculation from lit data
    int<lower=1> N_lit_ms; // Number of literature ms measurements
    int<lower=1> N_ms; // Number of ms measurements from our data

    // Observables for fitting
    vector<lower=0>[N_prot] protein_per_cell; 
    vector<lower=0>[N_prot] protein_per_cell_lambda; 
    vector<lower=1>[N_SA] SA;
    vector<lower=0>[N_SA] SA_lambda;
    vector<lower=0>[N_vol] cell_volume;
    vector<lower=0>[N_vol] cell_volume_lambda;
    vector<lower=0>[N_DNA] DNA_to_protein;

    // Observables for calculation
    vector<lower=0>[N_lit_ms] lit_ms_lam;
    vector<lower=0>[N_lit_ms] lit_phi_rib;
    vector<lower=0>[N_lit_ms] lit_phi_cyt;
    vector<lower=0>[N_lit_ms] lit_phi_mem;

    vector<lower=0>[N_ms] ms_lam;
    vector<lower=0>[N_ms] phi_rib;
    vector<lower=0>[N_ms] phi_cyt;
    vector<lower=0>[N_ms] phi_mem;
    vector<lower=0>[N_ms] phi_peri;
    vector<lower=0>[N_ms] vol;
    vector<lower=0>[N_ms] sa;
}

transformed data {
    // Log transforms
    vector[N_prot] log_prot_per_cell = log(protein_per_cell);
    vector[N_vol] log_volume = log(cell_volume);

}

parameters {
    // Linear regression on log protein per cell
    real beta0_prot; // Intercept
    real k_prot; // Slope
    real<lower=0> sigma_prot; // Homoskedastic error

    // Linear regression on surface area
    real beta0_sa; // Intercept
    real beta1_sa; // Slope
    real<lower=0> sigma_sa; // Homoskedastic error

    // Linear regression on log volume
    real beta0_vol; // intercept
    real beta1_vol; // slope
    real<lower=0> sigma_vol; // Homoskedastic error

    // Parameters for constant theta
    real<lower=0> theta_mu; // mean
    real<lower=0> theta_sigma; // standard deviation
}

model { 
    // Priors for total protein inference
    beta0_prot ~ normal(0, 1); 
    k_prot ~ normal(3, 2);
    sigma_prot ~ normal(0, 0.1);

    // Priors for surface area inference
    beta0_sa ~ normal(0, 3);
    beta1_sa ~ normal(0, 2); 
    sigma_sa ~ normal(0, 1);

    // Priors for volume inference
    beta0_vol ~ normal(0, 1);
    beta1_vol ~ normal(0, 1);
    sigma_vol ~ normal(0, 1);

    // Priors for DNA to protein ratio
    theta_mu ~ normal(0, 0.1);    
    theta_sigma ~ normal(0, 0.1);

    // Likelihoods 
    log_prot_per_cell ~ normal(beta0_prot + k_prot * protein_per_cell_lambda, sigma_prot);
    log_volume ~ normal(beta0_vol + beta1_vol * cell_volume_lambda, sigma_vol);
    SA ~ normal(beta0_sa + beta1_sa * SA_lambda, sigma_sa);
    DNA_to_protein ~ normal(theta_mu, theta_sigma);
}

generated quantities {
    // Calculate the quantities from the lit mass spec data
    vector<lower=0>[N_lit_ms] lit_ms_prot = exp(beta0_prot + k_prot * lit_ms_lam);
    vector<lower=0>[N_lit_ms] lit_ms_vol = exp(beta0_vol + beta1_vol * lit_ms_lam);
    vector<lower=0>[N_lit_ms] lit_ms_sa = beta0_sa + beta1_sa * lit_ms_lam;

    // Compute the lam scan properties
    vector<lower=0>[N_scan] scan_vol = exp(beta0_vol + beta1_vol * scan_lam);
    vector<lower=0>[N_scan] scan_sa = beta0_sa + beta1_sa * scan_lam;
    vector<lower=0>[N_scan] scan_prot = exp(beta0_prot + k_prot * scan_lam);
    vector<lower=0>[N_scan] scan_theta;
    for (i in 1:N_scan) {
        scan_theta[i] = theta_mu;
    }

    // Calculate quantities from our ms data
    vector<lower=0>[N_ms] ms_prot = exp(beta0_prot + k_prot * ms_lam);
    
    // Calculate the density and ratios for literature
    vector<lower=0>[N_lit_ms] lit_rho_cyt = lit_ms_prot .* (2.19 .* lit_phi_rib + lit_phi_cyt + theta_mu) ./  (lit_ms_vol - 0.025 .* lit_ms_sa);
    vector<lower=0>[N_lit_ms] lit_sigma_mem = lit_ms_prot .* lit_phi_mem ./ (2 * lit_ms_sa);
    vector<lower=0>[N_lit_ms] lit_density_ratio = lit_rho_cyt ./ lit_sigma_mem;

    // Calculate the density and ratios for our data
    vector<lower=0>[N_ms] rho_cyt = ms_prot .* (2.19 .* phi_rib + phi_cyt + theta_mu) ./  (vol - 0.025 .* sa);
    vector<lower=0>[N_ms] sigma_mem = ms_prot .* phi_mem ./ (2 * sa);
    vector<lower=0>[N_ms] density_ratio = rho_cyt ./ sigma_mem;

    // For our data, compute the predicted SAV given the density ratio
    vector<lower=0>[N_ms] predicted_SAV = (phi_mem .* density_ratio) ./ (2 .* (1 + 2.19 * phi_rib - phi_mem - phi_peri));

}