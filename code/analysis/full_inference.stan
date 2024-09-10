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

    // our mass spec measurements for for calculation 
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
    }

transformed data {
    // Log transforms of data with exponential fits
    vector<lower=1>[N_prot] log_prot_per_cell = log(prot_per_cell);
    vector[N_size] log_volume = log(volume);
    vector[N_size] log_sa = log(sa);
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
}

model {
   // Define the priors
   beta_0_prot ~ normal(0, 1);
   beta_1_prot ~ normal(0, 1);
   prot_sigma ~ normal(0, 1);
   beta_0_vol ~ normal(0, 1);
   beta_1_vol ~ normal(0, 1);
   vol_sigma ~ normal(0, 1);
   beta_0_sa ~ normal(0, 1);
   beta_1_sa ~ normal(0, 1);
   sa_sigma ~ normal(0, 1);

   // Define the likelihoods
   log_prot_per_cell ~ normal(beta_0_prot + beta_1_prot * prot_growth_rate_hr, prot_sigma);
   log_volume ~ normal(beta_0_vol + beta_1_vol * vol_growth_rate_hr, vol_sigma);
   log_sa ~ normal(beta_0_sa + beta_1_sa * sa_growth_rate_hr, sa_sigma);
}

generated quantities {
    // Define the quantities to calculate for the literature data
    vector[N_obs_lit] lit_rho_peri_ppc;
    vector[N_obs_lit] lit_rho_cyto_ppc;
    vector[N_obs_lit] lit_sigma_mem_ppc;
    vector[N_obs_lit] lit_prot_ppc;
    vector[N_obs_lit] lit_volume_ppc;
    vector[N_obs_lit] lit_sa_ppc;
    vector[N_obs_lit] lit_m_peri_ppc;
    vector[N_obs_lit] lit_kappa_ppc;
    real lit_kappa_ppc_mean;

    // Define the quantities to calculate for the observed data
    vector[N_obs] rho_peri_ppc;
    vector[N_obs] rho_cyto_ppc;
    vector[N_obs] sigma_mem_ppc;
    vector[N_obs] prot_ppc;
    vector[N_obs] m_peri_ppc;
    vector[N_obs] kappa_ppc;
    real kappa_ppc_mean;
    real weighted_kappa_ppc;

    // Define the quantities to calculate for the fits
    vector[N_pred] fit_prot;
    vector[N_pred] fit_volume;
    vector[N_pred] fit_sa;


    // Calculate the quantites for the literature data
    for (i in 1:N_obs_lit) {
        // Infer the necessary parameters
        lit_prot_ppc[i] = exp(normal_rng(beta_0_prot + beta_1_prot * lit_growth_rate_hr[i], prot_sigma));
        lit_volume_ppc[i] = exp(normal_rng(beta_0_vol + beta_1_vol * lit_growth_rate_hr[i], vol_sigma));
        lit_sa_ppc[i] = exp(normal_rng(beta_0_sa + beta_1_sa * lit_growth_rate_hr[i], sa_sigma));

        // Calculate the densities
        lit_rho_peri_ppc[i] = lit_phi_peri[i] * lit_prot_ppc[i] / (DELTA_PERI * lit_sa_ppc[i]);
        lit_sigma_mem_ppc[i] = lit_phi_mem[i] * lit_prot_ppc[i] / (2 * lit_sa_ppc[i]);
        lit_rho_cyto_ppc[i] = lit_prot_ppc[i] * (BETA_RIB * lit_phi_rib[i] + lit_phi_cyto[i]) / (lit_volume_ppc[i] - DELTA_PERI *lit_sa_ppc[i]);
        lit_kappa_ppc[i] = lit_rho_cyto_ppc[i] / lit_sigma_mem_ppc[i];
        lit_m_peri_ppc[i] = lit_phi_peri[i] * lit_prot_ppc[i];

    }
    lit_kappa_ppc_mean = mean(lit_kappa_ppc);
    
// Calculate the quantiteis for the observed data
    for (i in 1:N_obs) {
        prot_ppc[i] = exp(normal_rng(beta_0_prot + beta_1_prot * obs_growth_rate_hr[i], prot_sigma));
        rho_peri_ppc[i] = obs_phi_peri[i] * prot_ppc[i] / (DELTA_PERI * obs_sa[i]);
        rho_cyto_ppc[i] = prot_ppc[i] * (BETA_RIB * obs_phi_rib[i] + obs_phi_cyto[i]) / (obs_volume[i] - DELTA_PERI * obs_sa[i]);
        sigma_mem_ppc[i] = obs_phi_mem[i] * prot_ppc[i] / (2 * obs_sa[i]);
        m_peri_ppc[i] = obs_phi_peri[i] * prot_ppc[i];
        kappa_ppc[i] = rho_cyto_ppc[i] / sigma_mem_ppc[i];
    }
    kappa_ppc_mean = mean(kappa_ppc);

    weighted_kappa_ppc = ((kappa_ppc_mean/variance(kappa_ppc)) + (lit_kappa_ppc_mean/variance(lit_kappa_ppc)))/(1/variance(kappa_ppc) + 1/variance(lit_kappa_ppc));
    // Calculate the continuous PPCs for the fits
    for (i in 1:N_pred) {
        fit_prot[i] = exp(normal_rng(beta_0_prot + beta_1_prot * pred_growth_rate_hr[i], prot_sigma));
        fit_volume[i] = exp(normal_rng(beta_0_vol + beta_1_vol * pred_growth_rate_hr[i], vol_sigma));
        fit_sa[i] = exp(normal_rng(beta_0_sa + beta_1_sa * pred_growth_rate_hr[i], sa_sigma));
    }

}
