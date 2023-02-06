data {
    // Total condition dimensionality
    int<lower=1> J_size_cond; // Number of unique conditions with size measurements
    int<lower=1> J_growth_cond;  //  Number of unique growth conditions
    int<lower=1> J_growth_curves; // Number of unique growth curves
    int<lower=1> J_brad_cond;  // Number of unique conditions in bradford measurements
    int<lower=1> J_flow_cond; // Number of unique conditions for flow measurements 

    // Total measurement dimensionality
    int<lower=1> N_brad; // Number of bradford measurements
    int<lower=1> N_flow; // Number of flow measurements
    int<lower=1> N_growth; // Total number of growth measurements
    int<lower=1> N_size; // Total number of size measurements

    // Set linking id vectors 
    array[N_brad] int<lower=1, upper=J_size_cond> brad_link_idx;
    array[N_flow] int<lower=1, upper=J_size_cond> flow_link_idx;
    array[J_growth_curves] int<lower=1, upper=J_size_cond> growth_link_idx;
    array[N_size] int<lower=1, upper=J_size_cond> size_idx;


    // Set individual element dimensionality
    array[N_brad] int<lower=1, upper=J_brad_cond> brad_idx;
    array[N_flow] int<lower=1, upper=J_flow_cond> flow_idx;
    array[J_growth_curves] int<lower=1, upper=J_growth_cond> growth_cond_idx;
    array[N_growth] int<lower=1, upper=J_growth_curves> growth_curve_idx;

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
    //------------------------------------------------------------------------- 
    // Measured information
    vector<lower=0>[N_brad] prot_od; // OD595 per unit biomass measurements for periplasmic protein
    vector<lower=0>[N_brad] od_conv_factor; // Net conversion factor accounting for volume diffeerences

    // -------------------------------------------------------------------------
    // Literature Values for dry mass and cell counts
    // -------------------------------------------------------------------------
    int<lower=1> N_lit_meas; 
    vector<lower=1>[N_lit_meas] drymass;

    // -------------------------------------------------------------------------
    // Growth rate measurements
    // -------------------------------------------------------------------------
    vector<lower=0>[N_growth] growth_time; // Elapsed time for growth curves
    vector<lower=0>[N_growth] growth_od; // OD600nm measurements for growth curves

    // -------------------------------------------------------------------------
    // Flow cytometry measurements
    // -------------------------------------------------------------------------    
    // Measured information
    vector<lower=1>[N_flow] cells_per_biomass; // Measurements of  cells per mL of experimental culture

    // -------------------------------------------------------------------------
    // Size Data
    // -------------------------------------------------------------------------
    // Measured parameters
    vector<lower=0>[N_size] mean_width;
    vector<lower=0>[N_size] mean_length;
    vector<lower=0>[N_size] mean_vol;
    vector<lower=0>[N_size] mean_peri_vol;
    vector<lower=0>[N_size] mean_peri_vol_frac;
    vector<lower=0>[N_size] mean_sa;
    vector<lower=0>[N_size] mean_sav;

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

//     // -------------------------------------------------------------------------
//     // Linking Indices
//     // -------------------------------------------------------------------------
//     array[N_cell_meas] int<lower=1, upper=J_growth_cond> cells_per_biomass_growth_idx; 
//     array[J_prot_cond] int<lower=1, upper=J_growth_cond> prot_cond_map;
}

transformed data{
    // -------------------------------------------------------------------------
    // Protein Quantification Transformations 
    // -------------------------------------------------------------------------
    // Standardize input data 
    vector[N_brad] log_prot_od = log(prot_od);

    // -------------------------------------------------------------------------
    // Lit data transformations
    // -------------------------------------------------------------------------
    vector[N_lit_meas] drymass_tilde = (drymass - mean(drymass)) ./ sd(drymass);

    // -------------------------------------------------------------------------
    // Growth Rate Transformations
    // -------------------------------------------------------------------------
    // vector[N_growth]  log_growth_od = log(growth_od); // Compute the log density to permit linear regression
    
    // -------------------------------------------------------------------------
    // Flow Cytometry Transformations
    // -------------------------------------------------------------------------
    // vector[N_flow] log_billion_cells_per_biomass = log(1 / (cells_per_biomass ./ 1E9));

    // -------------------------------------------------------------------------
    // Size Measurement Transformations
    // -------------------------------------------------------------------------
    vector[N_size] mean_width_tilde = (mean_width - mean_width_mean[size_idx]) ./ mean_width_std[size_idx];
    vector[N_size] mean_length_tilde = (mean_length - mean_length_mean[size_idx]) ./ mean_length_std[size_idx];
    vector[N_size] mean_vol_tilde = (mean_vol - mean_vol_mean[size_idx]) ./ mean_vol_std[size_idx];
    vector[N_size] mean_peri_vol_tilde = (mean_peri_vol - mean_peri_vol_mean[size_idx]) ./ mean_peri_vol_std[size_idx];
    vector[N_size] mean_peri_vol_frac_tilde = (mean_peri_vol_frac - mean_peri_vol_frac_mean[size_idx]) ./ mean_peri_vol_frac_std[size_idx];
    vector[N_size] mean_sa_tilde = (mean_sa - mean_sa_mean[size_idx]) ./ mean_sa_std[size_idx];
    vector[N_size] mean_sav_tilde = (mean_sav - mean_sav_mean[size_idx]) ./ mean_sav_std[size_idx];    
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
    vector[J_brad_cond] log_prot_per_biomass_mu;
    vector<lower=0>[J_brad_cond] od595_per_biomass_sigma;

    // -------------------------------------------------------------------------
    // Lit data parameters
    // -------------------------------------------------------------------------
    real drymass_mu_tilde;
    real<lower=0> drymass_sigma_tilde;
   
    // // -------------------------------------------------------------------------
    // // Growth Rate Parameters
    // // -------------------------------------------------------------------------
    // // Hyper parameters 
    // vector<lower=0>[J_growth_cond] growth_mu;
    // real growth_tau;

    // // Level 1 hyper parameters
    // vector<lower=0>[J_growth_curves] growth_mu_1_tilde;

    // // Singular
    // real<lower=0> growth_sigma;
    // vector<lower=0>[J_growth_curves] growth_od_init;

    // -------------------------------------------------------------------------
    // Flow Cytometry Parameters
    // -------------------------------------------------------------------------
    // real biomass_per_volume;  
    // real biomass_per_vol_intercept;
    // real<lower=0> cells_per_biomass_sigma;

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
    // // -------------------------------------------------------------------------
    // // Protein Quantification Transformed Parameters
    // // -------------------------------------------------------------------------
    vector<lower=0>[J_brad_cond] prot_per_biomass_mu = exp(log_prot_per_biomass_mu);

    // -------------------------------------------------------------------------
    // Literature data transformed parameters
    // -------------------------------------------------------------------------
    // real<lower=0> drymass_mu = drymass_mu_tilde * sd(drymass) + mean(drymass);

    // -------------------------------------------------------------------------
    // Flow data transformed parameters
    // -------------------------------------------------------------------------
    // real  flow_slope = biomass_per_volume * drymass_mu;    

    // // -------------------------------------------------------------------------
    // // Growth Rate Transformed Parameters
    // // -------------------------------------------------------------------------
    // vector[J_growth_curves] growth_mu_1 = growth_mu[growth_cond_idx] + growth_tau * growth_mu_1_tilde;
    // vector[J_growth_curves] log_growth_od_init = log(growth_od_init);

    // // -------------------------------------------------------------------------
    // Size Measurement Transformed Parameters
    // -------------------------------------------------------------------------
    vector[J_size_cond] width_mu = mean_width_std .* width_mu_tilde + mean_width_mean;
    vector[J_size_cond] length_mu = mean_length_std .* length_mu_tilde + mean_length_mean;
    vector[J_size_cond] vol_mu = mean_vol_std .* vol_mu_tilde + mean_vol_mean;
    vector[J_size_cond] peri_vol_mu = mean_peri_vol_std .* peri_vol_mu_tilde + mean_peri_vol_mean;
    vector[J_size_cond] peri_vol_frac_mu = mean_peri_vol_frac_std .* peri_vol_frac_mu_tilde + mean_peri_vol_frac_mean;
    vector[J_size_cond] sa_mu = mean_sa_std .* sa_mu_tilde + mean_sa_mean;
    vector[J_size_cond] sav_mu = mean_sav_std .* sav_mu_tilde + mean_sav_mean;

    // -------------------------------------------------------------------------
    // Cells per Biomass Transformed Parameters
    // -------------------------------------------------------------------------
    // vector<lower=0>[J_growth_cond] vol_cells_per_biomass = 1E9 .* (beta_0_tilde - k_cells_per_biomass_tilde * vol_mu);



}

model {
    // -------------------------------------------------------------------------
    // Protein Quantification Model
    // -------------------------------------------------------------------------
    // Hyper priors
    cal_slope ~ normal(0, 0.1);
    cal_intercept ~ normal(0, 0.1);


    // Low-level priors
    cal_sigma ~ normal(0, 1);

    // Measurement priors
    log_prot_per_biomass_mu ~ std_normal();
    od595_per_biomass_sigma ~ std_normal();

    // Likelihoods
    cal_od ~ normal(cal_intercept + cal_slope .* concentration, cal_sigma);
    log_prot_od ~ cauchy(log(cal_intercept + cal_slope * prot_per_biomass_mu[brad_idx]./od_conv_factor), od595_per_biomass_sigma[brad_idx]);

    // // -------------------------------------------------------------------------
    // // Literature data model
    // // -------------------------------------------------------------------------
    // drymass_mu_tilde ~ std_normal();
    // drymass_sigma_tilde ~ std_normal();
    // drymass_tilde ~ normal(drymass_mu_tilde, drymass_sigma_tilde);

    // // -------------------------------------------------------------------------
    // // Growth Rate Model
    // // -------------------------------------------------------------------------
    // // Priors
    // growth_mu ~ gamma(2.85, 3.30);
    // growth_tau ~ std_normal();
    // growth_mu_1_tilde ~ std_normal();
    // growth_sigma ~ normal(0, 0.1); 
    // growth_od_init ~ normal(0, 0.1);

    // // Likelihood
    // log_growth_od ~ normal(log_growth_od_init[growth_curve_idx] + growth_mu_1[growth_curve_idx] .* growth_time, growth_sigma);

    // -------------------------------------------------------------------------
    // Flow cytometry model
    // -------------------------------------------------------------------------
    // Priors
    // biomass_per_volume ~ normal(0, 0.1);
    // cells_per_biomass_sigma ~ std_normal();

    // Likelihood
    // log_billion_cells_per_biomass ~ normal(log(flow_slope .* vol_mu[flow_link_idx]), cells_per_biomass_sigma);
    // log_billion_cells_per_biomass ~ cauchy(log(k_cells_per_biomass_tilde * vol_mu[cells_per_biomass_growth_idx]), cells_per_biomass_sigma);
 
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
    vector[N_brad] od595_per_biomass_rep;
    vector[N_cal_meas] od595_calib_rep; 
    //Posterior Predictive Checks
    for (i in 1:N_brad) {
        od595_per_biomass_rep[i] = exp(cauchy_rng(log(cal_intercept + cal_slope * prot_per_biomass_mu[brad_idx[i]] / od_conv_factor[i]), od595_per_biomass_sigma[brad_idx[i]])) * od_conv_factor[i];
    }
    for (i in 1:N_cal_meas) { 
        od595_calib_rep[i] = normal_rng(cal_intercept + cal_slope * concentration[i] , cal_sigma);
    } 

    // // -------------------------------------------------------------------------
    // // Lit data ppc
    // // -------------------------------------------------------------------------
    // vector[N_lit_meas] drymass_rep;
    // for (i in 1:N_lit_meas) {
    //     drymass_rep[i] = normal_rng(drymass_mu_tilde, drymass_sigma_tilde) * sd(drymass) + mean(drymass);
    // }

    // // -------------------------------------------------------------------------
    // // Growth Rate PPC
    // // -------------------------------------------------------------------------
    // vector[N_growth_meas] growth_od_rep;
    // for (i in 1:N_growth_meas) {
    //     growth_od_rep[i] = exp(normal_rng(log_growth_od_init[growth_curve_idx[i]] + growth_mu_1[growth_curve_idx[i]] .* growth_time[i], growth_sigma));
    // }

    // -------------------------------------------------------------------------
    // Cells Per Biomass PPC
    // -------------------------------------------------------------------------
    // vector[N_flow] cells_per_biomass_rep;
    // for (i in 1:N_flow) {
    //     cells_per_biomass_rep[i] = 1E9 * exp(1 / normal_rng(log(flow_slope * vol_mu[flow_link_idx[i]]), cells_per_biomass_sigma));
        
        // og(beta_0_tilde - k_cells_per_biomass_tilde * vol_mu[cells_per_biomass_growth_idx[i]]), cells_per_biomass_sigma));
    // }

    // -------------------------------------------------------------------------
    // Size PPC
    // -------------------------------------------------------------------------
    vector[N_size] width_rep;
    vector[N_size] length_rep;
    vector[N_size] vol_rep;
    vector[N_size] peri_vol_rep;
    vector[N_size] peri_vol_frac_rep;
    vector[N_size] sa_rep;
    vector[N_size] sav_rep;

    for (i in 1:N_size) {
        width_rep[i] = mean_width_mean[size_idx[i]] + mean_width_std[size_idx[i]] * normal_rng(width_mu_tilde[size_idx[i]], width_sigma[size_idx[i]]); 
        length_rep[i] = mean_length_mean[size_idx[i]] + mean_length_std[size_idx[i]] * normal_rng(length_mu_tilde[size_idx[i]], length_sigma[size_idx[i]]); 
        vol_rep[i] = mean_vol_mean[size_idx[i]] + mean_vol_std[size_idx[i]] * normal_rng(vol_mu_tilde[size_idx[i]], vol_sigma[size_idx[i]]); 
        peri_vol_rep[i] = mean_peri_vol_mean[size_idx[i]] + mean_peri_vol_std[size_idx[i]] * normal_rng(peri_vol_mu_tilde[size_idx[i]], peri_vol_sigma[size_idx[i]]); 
        peri_vol_frac_rep[i] = mean_peri_vol_frac_mean[size_idx[i]] + mean_peri_vol_frac_std[size_idx[i]] *normal_rng(peri_vol_frac_mu_tilde[size_idx[i]], peri_vol_frac_sigma[size_idx[i]]); 
        sa_rep[i] = mean_sa_mean[size_idx[i]] + mean_sa_std[size_idx[i]] * normal_rng(sa_mu_tilde[size_idx[i]], sa_sigma[size_idx[i]]); 
        sav_rep[i] = mean_sav_mean[size_idx[i]] + mean_sav_std[size_idx[i]] * normal_rng(sav_mu_tilde[size_idx[i]], sav_sigma[size_idx[i]]); 
    }

    // // -------------------------------------------------------------------------
    // // Compute properties
    // // -------------------------------------------------------------------------
    // vector[J_prot_cond] peri_density = prot_per_biomass_mu * 1E9 ./ (vol_cells_per_biomass[prot_cond_map] .* peri_vol_mu[prot_cond_map]); 
    // vector[J_prot_cond] peri_drymass_frac = prot_per_biomass_mu ./ drymass_mu;
    // vector[J_prot_cond] peri_protein_frac = prot_per_biomass_mu ./ (0.55 * drymass_mu);
}
