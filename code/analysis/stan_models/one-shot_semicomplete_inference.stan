data {

    // -------------------------------------------------------------------------
    // Bradford Calibration Curve
    // -------------------------------------------------------------------------
    // Dimensional information
    int<lower=1> J_cal_brep; // Number of biological replicates for Bradford calibration curve
    int<lower=1> N_cal_meas; // Number of Bradford calibration curve measurements
    array[N_cal_meas] int<lower=1, upper=J_cal_brep> cal_idx; // ID vector for calibration curve

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
    vector<lower=0>[J_prot_cond] mean_prot_od;
    vector<lower=0>[J_prot_cond] std_prot_od;

    // // -------------------------------------------------------------------------
    // // Flow cytometry measurements
    // // -------------------------------------------------------------------------
    // // Dimensional information
    // int<lower=1> J_flow_cond; // Number of conditions measured in flow cytometry data
    // int<lower=1> J_flow_trep; // Number of technical replicates
    // int<lower=1> J_flow_brep; // Number of biological replicates
    // int<lower=1> N_flow_meas; // Number of measurements
    // array[N_flow_meas] int<lower=1, upper=J_flow_trep> flow_trep_idx; // ID vector for technical replicates
    // array[J_flow_trep] int<lower=1, upper=J_flow_cond> flow_brep_idx; // ID vector for biological replicates
    // array[J_flow_brep] int<lower=1, upper=J_flow_cond> flow_cond_idx; // ID vector for conditions

    // // Measured information
    // vector[N_flow_meas]<lower=1> flow_meas; // Measurements of  cells per mL of experimental culture

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
    // Size Data
    // -------------------------------------------------------------------------
    // TODO

    // -------------------------------------------------------------------------
    // Condition linking vectors
    // -------------------------------------------------------------------------
    // TODO
}

transformed data{
    // -------------------------------------------------------------------------
    // Protein Quantification Transformations 
    // -------------------------------------------------------------------------
    // Standardize input data 
    vector[N_prot_meas] prot_od_tilde = (prot_od - mean_prot_od[prot_idx]) ./ (std_prot_od[prot_idx]);

    // -------------------------------------------------------------------------
    // Growth Rate Transformed Parameters
    // -------------------------------------------------------------------------
    vector[N_growth_meas]  log_growth_od = log(growth_od); // Compute the log density to permit linear regression
}

parameters {
    // -------------------------------------------------------------------------
    // Protein Quantification Parameters
    // -------------------------------------------------------------------------
    // Hyperparameters for calibration
    real<lower=0> cal_slope;
    real<lower=0> cal_intercept;
    real cal_tau;

    // Low-level parameter
    vector<lower=0>[J_cal_brep] cal_slope_1_tilde;
    vector<lower=0>[J_cal_brep] cal_intercept_1_tilde;
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
}

transformed parameters {
    // -------------------------------------------------------------------------
    // Protein Quantification Transformed Parameters
    // -------------------------------------------------------------------------
    vector<lower=0>[J_cal_brep] cal_slope_1 = cal_slope_1_tilde + cal_tau * cal_slope;
    vector<lower=0>[J_cal_brep] cal_intercept_1 = cal_intercept_1_tilde + cal_tau * cal_intercept;
    vector<lower=0>[J_prot_cond] od595_per_biomass_mu = (od595_per_biomass_mu_tilde .* std_prot_od) + mean_prot_od;
    vector<lower=0>[J_prot_cond] prot_per_biomass = (od595_per_biomass_mu - cal_intercept) ./ cal_slope;

    // -------------------------------------------------------------------------
    // Growth Rate Transformed Parameters
    // -------------------------------------------------------------------------
    vector[J_growth_brep] growth_mu_1 = growth_mu[growth_cond_idx] + growth_tau * growth_mu_1_tilde;
    vector[J_growth_brep] log_growth_od_init = log(growth_od_init);
}

model {
    // -------------------------------------------------------------------------
    // Protein Quantification Model
    // -------------------------------------------------------------------------
    // Hyper priors
    cal_slope ~ normal(0, 1);
    cal_intercept ~ normal(0, 1);
    cal_tau ~ std_normal();

    // Low-level priors
    cal_slope_1_tilde ~ std_normal();
    cal_intercept_1_tilde ~ std_normal();
    cal_sigma ~ normal(0, 1);

    // Measurement priors
    od595_per_biomass_mu_tilde ~ std_normal();
    od595_per_biomass_sigma ~ std_normal();//normal(0, 0.1);

    // Likelihoods
    cal_od ~ normal(cal_intercept_1[cal_idx] + cal_slope_1[cal_idx] .* concentration, cal_sigma);
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
        od595_calib_rep[i] = normal_rng(cal_intercept_1[cal_idx[i]] + cal_slope_1[cal_idx[i]] * concentration[i] , cal_sigma);
    } 

    // -------------------------------------------------------------------------
    // Growth Rate PPC
    // -------------------------------------------------------------------------
    vector[N_growth_meas] growth_od_rep;
    for (i in 1:N_growth_meas) {
        growth_od_rep[i] = exp(normal_rng(log_growth_od_init[growth_brep_idx[i]] + growth_mu_1[growth_brep_idx[i]] .* growth_time[i], growth_sigma));
    }

}
