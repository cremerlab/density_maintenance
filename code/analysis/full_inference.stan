data { 
    // Dimensionality
    int<lower=1> N_prot; // Number of protein per cell measurements data points
    int<lower=1> N_size; // Number of volume measurements data points (for lit MS measurements)
    int<lower=1> N_obs; // Number of our size and mass spec measurements 
    int<lower=1> N_obs_lit; // Number of literature size and mass spec measurements
    int<lower=1> N_pred; // Number of points for poseterior predictive check measurements

    // Literature data for estimating trends
    vector<lower=0>[N_prot] prot_growth_rate_hr; 
    vector<lower=0>[N_prot] prot_per_cell;
    vector<lower=0>[N_size] vol_growth_rate_hr;
    vector<lower=0>[N_size] volume;
    vector<lower=0>[N_size] sa_growth_rate_hr;
    vector<lower=0>[N_size] sa;

    // literature mass spec measurements for calculations    
    vector<lower=0>[N_obs_lit] lit_growth_rate_hr;
    vector<lower=0, upper=1>[N_obs_lit] lit_phi_peri;
    vector<lower=0, upper=1>[N_obs_lit] lit_phi_cyto;
    vector<lower=0, upper=1>[N_obs_lit] lit_phi_mem;
    vector<lower=0, upper=1>[N_obs_lit] lit_phi_rib;

    // Our mass spec measurements for for calculation and test of theory
    vector<lower=0>[N_obs] obs_growth_rate_hr;
    vector<lower=0>[N_obs] obs_sa;
    vector<lower=0>[N_obs] obs_volume;
    vector<lower=0, upper=1>[N_obs] obs_phi_peri;
    vector<lower=0, upper=1>[N_obs] obs_phi_cyto;
    vector<lower=0, upper=1>[N_obs] obs_phi_mem;
    vector<lower=0, upper=1>[N_obs] obs_phi_rib;

    // Constants
    real<lower=0> DELTA_PERI; // the periplasmic width
    real<lower=0> BETA_RIB; // Conversion factor from ribosome fraction to total RNA 
    vector<lower=0>[N_pred] pred_growth_rate_hr;
    vector<lower=0, upper=1>[N_pred] pred_phi_rib;
    }

transformed data {
    // Log transforms of data with exponential fits
    vector<lower=1>[N_prot] log_prot_per_cell = log(prot_per_cell);
    vector[N_size] log_volume = log(volume);
    vector[N_size] log_sa = log(sa);
    vector[N_obs] log_obs_phi_peri = log(obs_phi_peri);
}

parameters {   
    // Parameters for trend fits
    real beta_0_prot;
    real beta_1_prot;
    real<lower=0> prot_sigma;
    real beta_0_vol;
    real beta_1_vol;
    real<lower=0> vol_sigma;
    real beta_0_sa;
    real beta_1_sa;
    real<lower=0> sa_sigma;

    // Parameters for phi_rib dependence
    real<lower=0, upper=1> beta_0_phi_mem; // Intercept for membrane linear fit
    real beta_1_phi_mem; // Slope for membrane linear fit
    real<lower=0> phi_mem_sigma; // sigma for membrane linear fit
    real beta_0_phi_peri; // Intercept for periplasmic linear fit on log scale
    real beta_1_phi_peri; // Slope for periplasmic linear fit on log scale
    real<lower=0> phi_peri_sigma; // sigma for periplasmic linear fit on log scale
}

model {
   // Define the priors for trend fits
   beta_0_prot ~ normal(0, 1);
   beta_1_prot ~ normal(0, 1);
   prot_sigma ~ normal(0, 1);
   beta_0_vol ~ normal(0, 1);
   beta_1_vol ~ normal(0, 1);
   vol_sigma ~ normal(0, 1);
   beta_0_sa ~ normal(0, 1);
   beta_1_sa ~ normal(0, 1);
   sa_sigma ~ normal(0, 1);


   // Define the priors for phi_rib dependence. Quantiles set using distribution-explorer.github.io
   beta_0_phi_mem ~ beta(1.797, 17.74); // Quantiles [2.5, 97.5] = [0.01, 0.25]
   beta_1_phi_mem ~ normal(0, 1);
   phi_mem_sigma ~ normal(0, 0.1);
   beta_0_phi_peri ~ normal(0, 1);
   beta_1_phi_peri ~ normal(0, 10);
   phi_peri_sigma ~ normal(0, 0.1);

   // Define the likelihoods for size relations
   log_prot_per_cell ~ normal(beta_0_prot + beta_1_prot * prot_growth_rate_hr, prot_sigma);
   log_volume ~ normal(beta_0_vol + beta_1_vol * vol_growth_rate_hr, vol_sigma);
   log_sa ~ normal(beta_0_sa + beta_1_sa * sa_growth_rate_hr, sa_sigma);

   // Define the likelihoods for phi_rib dependence
   obs_phi_mem ~ normal(beta_0_phi_mem + beta_1_phi_mem * obs_phi_rib, phi_mem_sigma);
   log_obs_phi_peri ~ normal(beta_0_phi_peri + beta_1_phi_peri * obs_phi_rib, phi_peri_sigma);
}

generated quantities {
    // Define the quantities to calculate for the literature data
    vector[N_obs_lit] lit_rho_peri_ppc;
    vector[N_obs_lit] lit_rho_peri;
    vector[N_obs_lit] lit_rho_cyto_ppc;
    vector[N_obs_lit] lit_rho_cyto;
    vector[N_obs_lit] lit_sigma_mem_ppc;
    vector[N_obs_lit] lit_sigma_mem;
    vector[N_obs_lit] lit_prot_ppc;
    vector[N_obs_lit] lit_prot;
    vector[N_obs_lit] lit_volume_ppc;
    vector[N_obs_lit] lit_volume;
    vector[N_obs_lit] lit_sa_ppc;
    vector[N_obs_lit] lit_sa;
    vector[N_obs_lit] lit_m_peri_ppc;
    vector[N_obs_lit] lit_m_peri;
    vector[N_obs_lit] lit_kappa_ppc;
    vector[N_obs_lit] lit_kappa;
    real lit_kappa_ppc_mean;
    real lit_kappa_mean;

    // Define the quantities to calculate for the observed data
    vector[N_obs] rho_peri_ppc;
    vector[N_obs] rho_peri;
    vector[N_obs] rho_cyto_ppc;
    vector[N_obs] rho_cyto;
    vector[N_obs] sigma_mem_ppc;
    vector[N_obs] sigma_mem;
    vector[N_obs] prot_ppc;
    vector[N_obs] prot;
    vector[N_obs] m_peri_ppc;
    vector[N_obs] m_peri;
    vector[N_obs] kappa_ppc;
    vector[N_obs] kappa;
    vector[N_obs] pred_sav_ppc;
    vector[N_obs] pred_sav;
    real kappa_ppc_mean;
    real kappa_mean;
    real weighted_kappa_ppc;
    real weighted_kappa;

    // Define the quantities to calculate for the fits
    vector[N_pred] fit_prot;
    vector[N_pred] fit_volume;
    vector[N_pred] fit_sa;

    // Define quantities for predictions
    vector[N_pred] fit_phi_mem_ppc;
    vector[N_pred] fit_phi_mem;
    vector[N_pred] fit_phi_peri_ppc;
    vector[N_pred] fit_phi_peri;
    vector[N_pred] theory_sav_ppc;
    vector[N_pred] theory_sav;

    // Calculate the quantites for the literature data
    for (i in 1:N_obs_lit) {
        // Compute the posterior parameter distributions
        lit_prot[i] = exp(beta_0_prot + beta_1_prot * lit_growth_rate_hr[i]);
        lit_volume[i] = exp(beta_0_vol + beta_1_vol * lit_growth_rate_hr[i]);
        lit_sa[i] = exp(beta_0_sa + beta_1_sa * lit_growth_rate_hr[i]);

        // Compute the posterior predictive check
        lit_prot_ppc[i] = exp(normal_rng(log(lit_prot[i]), prot_sigma));
        lit_volume_ppc[i] = exp(normal_rng(log(lit_volume[i]), vol_sigma));
        lit_sa_ppc[i] = exp(normal_rng(log(lit_sa[i]), sa_sigma));

        // Calculate the true posterior densities and the posterior predictive distributions of each MS point
        lit_rho_peri[i] = lit_phi_peri[i] * lit_prot[i] / (DELTA_PERI * lit_sa[i]);
        lit_rho_peri_ppc[i] = lit_phi_peri[i] * lit_prot_ppc[i] / (DELTA_PERI * lit_sa_ppc[i]);
        lit_sigma_mem[i] = lit_phi_mem[i] * lit_prot[i] / (2 * lit_sa[i]);
        lit_sigma_mem_ppc[i] = lit_phi_mem[i] * lit_prot_ppc[i] / (2 * lit_sa_ppc[i]);
        lit_rho_cyto[i] = lit_prot[i] * (BETA_RIB * lit_phi_rib[i] + lit_phi_cyto[i]) / (lit_volume[i] - DELTA_PERI * lit_sa[i]);
        lit_rho_cyto_ppc[i] = lit_prot_ppc[i] * (BETA_RIB * lit_phi_rib[i] + lit_phi_cyto[i]) / (lit_volume_ppc[i] - DELTA_PERI *lit_sa_ppc[i]);
        lit_kappa[i] = lit_rho_cyto[i] / lit_sigma_mem[i];
        lit_kappa_ppc[i] = lit_rho_cyto_ppc[i] / lit_sigma_mem_ppc[i];
        lit_m_peri[i] = lit_phi_peri[i] * lit_prot[i];
        lit_m_peri_ppc[i] = lit_phi_peri[i] * lit_prot_ppc[i];

    }
    lit_kappa_ppc_mean = mean(lit_kappa_ppc);
    lit_kappa_mean = mean(lit_kappa);

    // Calculate the quantities for the observed data
    for (i in 1:N_obs) {
        prot[i] = exp(beta_0_prot + beta_1_prot * obs_growth_rate_hr[i]);
        prot_ppc[i] = exp(normal_rng(log(prot[i]), prot_sigma));
        rho_peri_ppc[i] = obs_phi_peri[i] * prot_ppc[i] / (DELTA_PERI * obs_sa[i]);
        rho_peri[i] = obs_phi_peri[i] * prot[i] / (DELTA_PERI * obs_sa[i]);
        rho_cyto_ppc[i] = prot_ppc[i] * (BETA_RIB * obs_phi_rib[i] + obs_phi_cyto[i]) / (obs_volume[i] - DELTA_PERI * obs_sa[i]);
        rho_cyto[i] = prot[i] * (BETA_RIB * obs_phi_rib[i] + obs_phi_cyto[i]) / (obs_volume[i] - DELTA_PERI * obs_sa[i]);
        sigma_mem_ppc[i] = obs_phi_mem[i] * prot_ppc[i] / (2 * obs_sa[i]);
        sigma_mem[i] = obs_phi_mem[i] * prot[i] / (2 * obs_sa[i]);
        m_peri_ppc[i] = obs_phi_peri[i] * prot_ppc[i];
        m_peri[i] = obs_phi_peri[i] * prot[i];
        kappa_ppc[i] = rho_cyto_ppc[i] / sigma_mem_ppc[i];
        kappa[i] = rho_cyto[i] / sigma_mem[i];
    }

    kappa_ppc_mean = mean(kappa_ppc);
    kappa_mean = mean(kappa);

    weighted_kappa_ppc = ((kappa_ppc_mean/variance(kappa_ppc)) + (lit_kappa_ppc_mean/variance(lit_kappa_ppc)))/(1/variance(kappa_ppc) + 1/variance(lit_kappa_ppc));
    weighted_kappa = ((kappa_mean/variance(kappa)) + (lit_kappa_mean/variance(lit_kappa)))/(1/variance(kappa) + 1/variance(lit_kappa));

    // Given the estimate for kappa, calculate the predicted SAV.
    for (i in 1:N_obs) {
        pred_sav_ppc[i] = kappa_ppc_mean * obs_phi_mem[i] / (2 * (1 + BETA_RIB * obs_phi_rib[i] - obs_phi_mem[i] - obs_phi_peri[i]));
        pred_sav[i] = kappa_mean * obs_phi_mem[i] / (2 * (1 + BETA_RIB * obs_phi_rib[i] - obs_phi_mem[i] - obs_phi_peri[i]));
    }

    // Calculate the continuous PPCs for the fits and for the predictions given our MS data
    for (i in 1:N_pred) {
        fit_prot[i] = exp(normal_rng(beta_0_prot + beta_1_prot * pred_growth_rate_hr[i], prot_sigma));
        fit_volume[i] = exp(normal_rng(beta_0_vol + beta_1_vol * pred_growth_rate_hr[i], vol_sigma));
        fit_sa[i] = exp(normal_rng(beta_0_sa + beta_1_sa * pred_growth_rate_hr[i], sa_sigma));
        fit_phi_mem[i] = beta_0_phi_mem + beta_1_phi_mem * pred_phi_rib[i];
        fit_phi_mem_ppc[i] = normal_rng(fit_phi_mem[i], phi_mem_sigma);
        fit_phi_peri[i] = exp(normal_rng(beta_0_phi_peri + beta_1_phi_peri * pred_phi_rib[i], phi_peri_sigma));
        fit_phi_peri_ppc[i] = exp(normal_rng(log(fit_phi_peri[i]), phi_peri_sigma));
        theory_sav_ppc[i] = fit_phi_mem_ppc[i] * kappa_ppc_mean / (2 * (1 + BETA_RIB * pred_phi_rib[i] - fit_phi_mem_ppc[i] - fit_phi_peri_ppc[i]));
        theory_sav[i] = fit_phi_mem[i] * kappa_mean / (2 * (1 + BETA_RIB * pred_phi_rib[i] - fit_phi_mem[i] - fit_phi_peri[i]));
    }
}
