data {
    // Dimensional Parameters
    int<lower=1> N_obs; // Number of measurements
    int<lower=1> N_fit; // Number of points over which to compute the ppc.    

    // Observed data
    vector<lower=0>[N_obs] obs_phi_rib; // Ribosomal allocation
    vector<lower=0>[N_obs] obs_phi_mem; // Membrane allocation
    vector<lower=0>[N_obs] obs_phi_peri; // Periplasmic allocation
    vector<lower=0>[N_obs] obs_sav; // Surface-to-volume ratio

    // Input data overwhich to draw the ppcs.
    vector<lower=0>[N_fit] phi_rib_range;

    // Input constants
    real<lower=0> BETA;
}

parameters {
    // Theory parameters
    real<lower=0> kappa; // Density ratio
    real<lower=0> sigma; // Homoskedastic error for theory

    // Allocation trend parameters
    real<lower=0> beta_0_phi_mem; // Intercept parameter for membrane allocation trend
    real beta_1_phi_mem; // Slope paramter for membrane allocation trend
    real<lower=0> sigma_phi_mem; // Homoskedastic error for membrane allocation trend
    real<lower=0> beta_0_phi_peri; // Intercept paramter for periplasmic allocation trend
    real beta_1_phi_peri; // Scale parameter for periplasmic allocation trend
    real<lower=0> sigma_phi_peri; // Homoskedastic error for periplasmic allocation trend
}

model {

    // Define the theoretical prediction
    vector[N_obs] mu = kappa .* obs_phi_mem ./ (2 .* (1 + BETA * obs_phi_rib - obs_phi_mem - obs_phi_peri));

    // Set priors. Parameter bounds chosen using
    //https://distribution-explorer.github.io/ with quantile setter mode. 
    // Percentile bounds are set as a comment in [2.5%, 97.5%]
    kappa ~ inv_gamma(2.863, 140.1); // [20 inv µm, 250 inv µm]
    sigma ~ normal(0, 1);
    beta_0_phi_mem ~ beta(1.262, 5.967); // [0.01, 0.5]
    beta_1_phi_mem ~ normal(0, 1);
    sigma_phi_mem ~ normal(0, 0.1);
    beta_0_phi_peri ~ beta(1.262, 5.967); // [0.01, 0.5]
    beta_1_phi_peri ~ normal(0, 10);
    sigma_phi_peri ~ normal(0, 0.1);

    // Set likelihoods
    obs_phi_mem ~ normal(beta_0_phi_mem + beta_1_phi_mem * obs_phi_rib, sigma_phi_mem);
    obs_phi_peri ~ normal(beta_0_phi_peri * exp(beta_1_phi_peri * obs_phi_rib), sigma_phi_peri);
    obs_sav ~ normal(mu, sigma);

}

generated quantities {
    vector[N_fit] phi_mem_mu;
    vector[N_fit] phi_mem_ppc;
    vector[N_fit] phi_peri_mu;
    vector[N_fit] phi_peri_ppc;
    vector[N_fit] theory_mu;
    vector[N_fit] theory_ppc;

    for (i in 1:N_fit) {
        phi_mem_mu[i] = beta_0_phi_mem + beta_1_phi_mem * phi_rib_range[i];
        phi_peri_mu[i] = beta_0_phi_peri * exp(beta_1_phi_peri * phi_rib_range[i]);
        phi_mem_ppc[i] = normal_rng(phi_mem_mu[i], sigma_phi_mem);
        phi_peri_ppc[i] = normal_rng(phi_peri_mu[i], sigma_phi_peri);
        theory_mu[i] = kappa * phi_mem_mu[i] / (2 * (1 + BETA * phi_rib_range[i] - phi_mem_mu[i] - phi_peri_mu[i]));
        theory_ppc[i] = normal_rng(theory_mu[i], sigma);
    }

}