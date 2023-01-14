data {

    // -------------------------------------------------------------------------
    // Bradford Calibration Curve
    // -------------------------------------------------------------------------
    // Dimensional information
    int<lower=1> N_cal_meas; // Number of Bradford calibration curve measurements

    // Measured information
    vector<lower=0>[N_cal_meas] concentration; // Protein concentrations for cal curve
    vector<lower=0>[N_cal_meas] cal_od; // OD595nm measurements for calibration curve

    // -------------------------------------------------------------------------
    // Bradford Protein Measurements
    //-------------------------------------------------------------------------- 

    // Dimensional information
    int<lower=1> J_prot_cond; // Number of conditions for protein measurements
    int<lower=1> N_prot_meas; // Number of measurements of periplasmic protein
    array[N_prot_meas] int<lower=1, upper=J_prot_cond> prot_idx; // ID vector for bradford measurements

    // Measured information
    vector<lower=0>[N_prot_meas] prot_od; // OD595 per unit biomass measurements for periplasmic protein

    // Measured information for standardization
    vector<lower=0>[J_prot_cond] mean_prot_od;
    vector<lower=0>[J_prot_cond] std_prot_od;

    // -------------------------------------------------------------------------
    // Growth rate measurements
    // -------------------------------------------------------------------------
    // Dimensional information
    int<lower=1> J_growth_cond; // Number of conditions for growth rate measurements
    int<lower=1> J_growth_brep; // Number of biological replicates of growth measurements
    int<lower=1> N_growth_meas; // Total number of OD measurements for all replicates and conditions
    array[J_growth_brep] int<lower=1, upper=J_growth_cond> growth_cond_idx; // ID vector for mapping
    array[N_growth_meas] int<lower=1, upper=J_growth_brep> growth_brep_idx; // ID vector for growth biological replicates

    // Measured information
    vector<lower=0>[N_growth_meas] growth_time; // Elapsed time for growth curves
    vector<lower=0>[N_growth_meas] growth_od; // OD600nm measurements for growth curves

    // -------------------------------------------------------------------------
    // Flow cytometry measurements
    // -------------------------------------------------------------------------
    // Dimensional information
    int<lower=1> N_cell_meas; // Total number of flow cytometry measurements
    
    // Measured information
    vector<lower=1>[N_cell_meas] cells_per_biomass; // Measurements of  cells per mL of experimental culture

    // -------------------------------------------------------------------------
    // Size Data
    // -------------------------------------------------------------------------
    // Dimensional Parameters
    int<lower=1> J_size_cond; // Number of conditions in which size was measured.
    int<lower=1> N_size_meas; // Total number of size measurements
    array[N_size_meas] int<lower=1, upper=J_size_cond> size_idx; //ID vector for sizes.

    // Measured parameters
    vector<lower=0>[N_size_meas] mean_width;
    vector<lower=0>[N_size_meas] mean_length;
    vector<lower=0>[N_size_meas] mean_vol;
    vector<lower=0>[N_size_meas] mean_peri_vol;
    vector<lower=0>[N_size_meas] mean_peri_vol_frac;
    vector<lower=0>[N_size_meas] mean_sa;
    vector<lower=0>[N_size_meas] mean_sav;

    // Parameters for standardization
    vector<lower=0>[J_size_cond] mean_width_mean;
    vector<lower=0>[J_size_cond] mean_width_std;
    vector<lower=0>[J_size_cond] mean_length_mean;
    vector<lower=0>[J_size_cond] mean_length_std;
    vector<lower=0>[J_size_cond] mean_vol_mean;
    vector<lower=0>[J_size_cond] mean_vol_std;
    vector<lower=0>[J_size_cond] mean_peri_vol_mean;
    vector<lower=0>[J_size_cond] mean_peri_vol_std;
    vector<lower=0>[J_size_cond] mean_peri_vol_frac_mean;
    vector<lower=0>[J_size_cond] mean_peri_vol_frac_std;
    vector<lower=0>[J_size_cond] mean_sa_mean;
    vector<lower=0>[J_size_cond] mean_sa_std;
    vector<lower=0>[J_size_cond] mean_sav_mean;
    vector<lower=0>[J_size_cond] mean_sav_std;

    // -------------------------------------------------------------------------
    // Linking Indices
    // -------------------------------------------------------------------------
    array[N_cell_meas] int<lower=1, upper=J_growth_cond> cells_per_biomass_growth_idx; 
    array[J_prot_cond] int<lower=1, upper=J_growth_cond> prot_cond_map;



}

transformed data{
    // -------------------------------------------------------------------------
    // Protein Quantification Transformations 
    // -------------------------------------------------------------------------
    // Standardize input data 
    vector[N_prot_meas] prot_od_tilde = (prot_od - mean_prot_od[prot_idx]) ./ (std_prot_od[prot_idx]);

    // -------------------------------------------------------------------------
    // Growth Rate Transformations
    // -------------------------------------------------------------------------
    vector[N_growth_meas]  log_growth_od = log(growth_od); // Compute the log density to permit linear regression
    
    // -------------------------------------------------------------------------
    // Flow Cytometry Transformations
    // -------------------------------------------------------------------------
    vector[N_cell_meas] billion_cells_per_biomass = cells_per_biomass ./ 1E9;

    // -------------------------------------------------------------------------
    // Size Measurement Transformations
    // -------------------------------------------------------------------------
    vector[N_size_meas] mean_width_tilde = (mean_width - mean_width_mean[size_idx]) ./ mean_width_std[size_idx];
    vector[N_size_meas] mean_length_tilde = (mean_length - mean_length_mean[size_idx]) ./ mean_length_std[size_idx];
    vector[N_size_meas] mean_vol_tilde = (mean_vol - mean_vol_mean[size_idx]) ./ mean_vol_std[size_idx];
    vector[N_size_meas] mean_peri_vol_tilde = (mean_peri_vol - mean_peri_vol_mean[size_idx]) ./ mean_peri_vol_std[size_idx];
    vector[N_size_meas] mean_peri_vol_frac_tilde = (mean_peri_vol_frac - mean_peri_vol_frac_mean[size_idx]) ./ mean_peri_vol_frac_std[size_idx];
    vector[N_size_meas] mean_sa_tilde = (mean_sa - mean_sa_mean[size_idx]) ./ mean_sa_std[size_idx];
    vector[N_size_meas] mean_sav_tilde = (mean_sav - mean_sav_mean[size_idx]) ./ mean_sav_std[size_idx];    
}


parameters {
    // -------------------------------------------------------------------------
    // Protein Quantification Parameters
    // -------------------------------------------------------------------------
    // Hyperparameters for calibration
    real<lower=0> cal_slope;
    real<lower=0> cal_intercept;
    real<lower=0> cal_sigma;

    // Parameters for protein measurements 
    vector[J_prot_cond] od595_per_biomass_mu_tilde;
    vector<lower=0>[J_prot_cond] od595_per_biomass_sigma;

    // -------------------------------------------------------------------------
    // Growth Rate Parameters
    // -------------------------------------------------------------------------
    // Hyper parameters 
    vector<lower=0>[J_growth_cond] growth_mu;
    real growth_tau;

    // Level 1 hyper parameters
    vector<lower=0>[J_growth_brep] growth_mu_1_tilde;

    // Singular
    real<lower=0> growth_sigma;
    vector<lower=0>[J_growth_brep] growth_od_init;

    // -------------------------------------------------------------------------
    // Flow Cytometry Parameters
    // -------------------------------------------------------------------------
    real<lower=0> k_cells_per_biomass_tilde;  
    real<lower=0> beta_0_tilde;
    real<lower=0> cells_per_biomass_sigma;

    // -------------------------------------------------------------------------
    // Size parameters
    // -------------------------------------------------------------------------
    vector[J_size_cond] width_mu_tilde;
    vector<lower=0>[J_size_cond] width_sigma;
    vector[J_size_cond] length_mu_tilde;
    vector<lower=0>[J_size_cond] length_sigma;
    vector[J_size_cond] vol_mu_tilde;
    vector<lower=0>[J_size_cond] vol_sigma;
    vector[J_size_cond] peri_vol_mu_tilde;
    vector<lower=0>[J_size_cond] peri_vol_sigma;
    vector[J_size_cond] peri_vol_frac_mu_tilde;
    vector<lower=0>[J_size_cond] peri_vol_frac_sigma;
    vector[J_size_cond] sa_mu_tilde;
    vector<lower=0>[J_size_cond] sa_sigma;
    vector[J_size_cond] sav_mu_tilde;
    vector<lower=0>[J_size_cond] sav_sigma;

}

transformed parameters {
    // -------------------------------------------------------------------------
    // Protein Quantification Transformed Parameters
    // -------------------------------------------------------------------------
    vector<lower=0>[J_prot_cond] od595_per_biomass_mu = (od595_per_biomass_mu_tilde .* std_prot_od) + mean_prot_od;
    vector<lower=0>[J_prot_cond] prot_per_biomass = (od595_per_biomass_mu - cal_intercept) ./ cal_slope;

    // -------------------------------------------------------------------------
    // Growth Rate Transformed Parameters
    // -------------------------------------------------------------------------
    vector[J_growth_brep] growth_mu_1 = growth_mu[growth_cond_idx] + growth_tau * growth_mu_1_tilde;
    vector[J_growth_brep] log_growth_od_init = log(growth_od_init);

    // -------------------------------------------------------------------------
    // Cells per Biomass Transformed Parameters
    // -------------------------------------------------------------------------
    vector<lower=0>[J_growth_cond] growth_cells_per_biomass = 1E9 .* (beta_0_tilde * exp(-k_cells_per_biomass_tilde * growth_mu));

    // -------------------------------------------------------------------------
    // Size Measurement Transformed Parameters
    // -------------------------------------------------------------------------
    vector[J_size_cond] width_mu = mean_width_std .* width_mu_tilde + mean_width_mean;
    vector[J_size_cond] length_mu = mean_length_std .* length_mu_tilde + mean_length_mean;
    vector[J_size_cond] vol_mu = mean_vol_std .* vol_mu_tilde + mean_vol_mean;
    vector[J_size_cond] peri_vol_mu = mean_peri_vol_std .* peri_vol_mu_tilde + mean_peri_vol_mean;
    vector[J_size_cond] peri_vol_frac_mu = mean_peri_vol_frac_std .* peri_vol_frac_mu_tilde + mean_peri_vol_frac_mean;
    vector[J_size_cond] sa_mu = mean_sa_std .* sa_mu_tilde + mean_sa_mean;
    vector[J_size_cond] sav_mu = mean_sav_std .* sav_mu_tilde + mean_sav_mean;

}

model {
    // -------------------------------------------------------------------------
    // Protein Quantification Model
    // -------------------------------------------------------------------------
    // Hyper priors
    cal_slope ~ normal(0, 0.1);
    cal_intercept ~ normal(0, 0.1);

    // Low-level priors
    cal_sigma ~ normal(0, 0.1);

    // Measurement priors
    od595_per_biomass_mu_tilde ~ std_normal();
    od595_per_biomass_sigma ~ std_normal();//normal(0, 0.1);

    // Likelihoods
    cal_od ~ normal(cal_intercept + cal_slope .* concentration, cal_sigma);
    prot_od_tilde ~ cauchy(od595_per_biomass_mu_tilde[prot_idx], od595_per_biomass_sigma[prot_idx]);

    // -------------------------------------------------------------------------
    // Growth Rate Model
    // -------------------------------------------------------------------------
    // Priors
    growth_mu ~ gamma(2.85, 3.30);
    growth_tau ~ std_normal();
    growth_mu_1_tilde ~ std_normal();
    growth_sigma ~ normal(0, 0.1); 
    growth_od_init ~ normal(0, 0.1);

    // Likelihood
    log_growth_od ~ normal(log_growth_od_init[growth_brep_idx] + growth_mu_1[growth_brep_idx] .* growth_time, growth_sigma);

    // -------------------------------------------------------------------------
    // Flow cytometry model
    // -------------------------------------------------------------------------
    // Priors
    k_cells_per_biomass_tilde ~ std_normal();
    beta_0_tilde ~ std_normal();
    cells_per_biomass_sigma ~ std_normal();

    // Likelihood
    billion_cells_per_biomass ~ normal(beta_0_tilde * exp(-k_cells_per_biomass_tilde * growth_mu[cells_per_biomass_growth_idx]), cells_per_biomass_sigma);
 
    // -------------------------------------------------------------------------
    // Size model 
    // -------------------------------------------------------------------------
    // Size priors
    width_mu_tilde ~ std_normal();
    width_sigma ~ normal(0, 0.1);
    length_mu_tilde ~ std_normal();
    length_sigma ~ normal(0, 0.1);
    vol_mu_tilde ~ std_normal();
    vol_sigma ~ normal(0, 0.1);
    peri_vol_mu_tilde ~ std_normal();
    peri_vol_sigma ~ normal(0, 0.1);
    peri_vol_frac_mu_tilde ~ std_normal();
    peri_vol_frac_sigma ~ normal(0, 0.1);
    sa_mu_tilde ~ std_normal();
    sa_sigma ~ normal(0, 0.1);
    sav_mu_tilde ~ std_normal();
    sav_sigma ~ normal(0, 0.1);

    // Likelihoods
    mean_width_tilde ~ normal(width_mu_tilde[size_idx], width_sigma[size_idx]);
    mean_length_tilde ~ normal(length_mu_tilde[size_idx], length_sigma[size_idx]);
    mean_vol_tilde  ~ normal(vol_mu_tilde[size_idx],  vol_sigma[size_idx]);
    mean_peri_vol_tilde ~ normal(peri_vol_mu_tilde[size_idx], peri_vol_sigma[size_idx]);
    mean_peri_vol_frac_tilde ~ normal(peri_vol_frac_mu_tilde[size_idx], peri_vol_frac_sigma[size_idx]);
    mean_sa_tilde ~ normal(sa_mu_tilde[size_idx], sa_sigma[size_idx]);
    mean_sav_tilde ~ normal(sav_mu_tilde[size_idx], sav_sigma[size_idx]);
}

generated quantities { 
    // -------------------------------------------------------------------------
    // Protein Quantification PPC
    // -------------------------------------------------------------------------
    vector<lower=0>[N_prot_meas] od595_per_biomass_rep;
    vector<lower=0>[N_cal_meas] od595_calib_rep; 
    //Posterior Predictive Checks
    for (i in 1:N_prot_meas) {
        od595_per_biomass_rep[i] = mean_prot_od[prot_idx[i]]  + (std_prot_od[prot_idx[i]] * normal_rng(od595_per_biomass_mu_tilde[prot_idx[i]], od595_per_biomass_sigma[prot_idx[i]]));
    }
    for (i in 1:N_cal_meas) { 
        od595_calib_rep[i] = normal_rng(cal_intercept + cal_slope * concentration[i] , cal_sigma);
    } 

    // -------------------------------------------------------------------------
    // Growth Rate PPC
    // -------------------------------------------------------------------------
    vector[N_growth_meas] growth_od_rep;
    for (i in 1:N_growth_meas) {
        growth_od_rep[i] = exp(normal_rng(log_growth_od_init[growth_brep_idx[i]] + growth_mu_1[growth_brep_idx[i]] .* growth_time[i], growth_sigma));
    }

    // -------------------------------------------------------------------------
    // Cells Per Biomass PPC
    // -------------------------------------------------------------------------
    vector[N_cell_meas] cells_per_biomass_rep;
    for (i in 1:N_cell_meas) {
        cells_per_biomass_rep[i] = 1E9 * normal_rng(beta_0_tilde * exp(-k_cells_per_biomass_tilde * growth_mu[cells_per_biomass_growth_idx[i]]), cells_per_biomass_sigma);
    }

    // -------------------------------------------------------------------------
    // Size PPC
    // -------------------------------------------------------------------------
    vector[N_size_meas] width_rep;
    vector[N_size_meas] length_rep;
    vector[N_size_meas] vol_rep;
    vector[N_size_meas] peri_vol_rep;
    vector[N_size_meas] peri_vol_frac_rep;
    vector[N_size_meas] sa_rep;
    vector[N_size_meas] sav_rep;

    for (i in 1:N_size_meas) {
        width_rep[i] = mean_width_mean[size_idx[i]] + mean_width_std[size_idx[i]] * normal_rng(width_mu_tilde[size_idx[i]], width_sigma[size_idx[i]]); 
        length_rep[i] = mean_length_mean[size_idx[i]] + mean_length_std[size_idx[i]] * normal_rng(length_mu_tilde[size_idx[i]], length_sigma[size_idx[i]]); 
        vol_rep[i] = mean_vol_mean[size_idx[i]] + mean_vol_std[size_idx[i]] * normal_rng(vol_mu_tilde[size_idx[i]], vol_sigma[size_idx[i]]); 
        peri_vol_rep[i] = mean_peri_vol_mean[size_idx[i]] + mean_peri_vol_std[size_idx[i]] * normal_rng(peri_vol_mu_tilde[size_idx[i]], peri_vol_sigma[size_idx[i]]); 
        peri_vol_frac_rep[i] = mean_peri_vol_frac_mean[size_idx[i]] + mean_peri_vol_frac_std[size_idx[i]] *normal_rng(peri_vol_frac_mu_tilde[size_idx[i]], peri_vol_frac_sigma[size_idx[i]]); 
        sa_rep[i] = mean_sa_mean[size_idx[i]] + mean_sa_std[size_idx[i]] * normal_rng(sa_mu_tilde[size_idx[i]], sa_sigma[size_idx[i]]); 
        sav_rep[i] = mean_sav_mean[size_idx[i]] + mean_sav_std[size_idx[i]] * normal_rng(sav_mu_tilde[size_idx[i]], sav_sigma[size_idx[i]]); 
    }

    // -------------------------------------------------------------------------
    // Compute properties
    // -------------------------------------------------------------------------
    vector[J_prot_cond] n_cells = 1E9 .* beta_0_tilde .* exp(-k_cells_per_biomass_tilde * growth_mu[prot_cond_map]); 
    vector[J_prot_cond] peri_density = prot_per_biomass * 1E9 ./ (n_cells .* peri_vol_mu[prot_cond_map]); 
    vector[J_prot_cond] peri_mass_frac = prot_per_biomass ./ (n_cells .* vol_mu[prot_cond_map] * 500);
 
}
