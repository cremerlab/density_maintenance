data { 
    int<lower=1> J_cond; // Number of unique carbon sources
    int<lower=1> J_calib; // Number of calibration curve replicates
    int<lower=1> N_meas; // Total number of measurements 
    int<lower=1> N_calib; // Total number of calibration measurements
    array[N_meas] int<lower=1, upper=J_cond> meas_idx; // ID vector mapping protein measurements to carbon source 
    array[N_calib] int<lower=1, upper=J_calib> calib_idx; // ID vector for mapping calibration curves

    // Measured data for calibration
    vector<lower=0>[N_calib] concentration;
    vector<lower=0>[N_calib] od_595nm_calib;

    // Measured data for quantification
    vector<lower=0>[N_meas] od_595nm_meas;
    vector<lower=0>[N_meas] od_600nm_meas;
    vector<lower=0>[N_meas] dil_factor;
    vector<lower=0>[N_meas] ext_volume;
    vector<lower=0>[N_meas] cult_volume;
    }

transformed data {
    vector<lower=0>[N_meas] od_per_biomass = od_595nm_meas .* dil_factor .* ext_volume ./ (cult_volume .* od_600nm_meas);
}

parameters {
    // Hyperparameters for calibration
    real slope_mu;
    real<lower=0> slope_sigma;
    real intercept_mu;
    real<lower=0> intercept_sigma;

    // Low-level parameter
    vector[J_calib] slope_1;
    vector[J_calib] intercept_1;
    real<lower=0> calib_sigma;

    // Parameters for protein measurements 
    vector<lower=0>[J_cond] od_per_biomass_mu;
    vector<lower=0>[J_cond] od_per_biomass_sigma;
}

transformed parameters {
    vector<lower=0>[J_cond] prot_per_biomass = (od_per_biomass_mu - intercept_mu) ./ slope_mu;
}

model {
    // Hyper priors
    slope_mu ~ std_normal();
    slope_sigma ~ std_normal(); 
    intercept_mu ~ std_normal();
    intercept_sigma ~ std_normal();

    // Low-level priors
    slope_1 ~ normal(slope_mu, slope_sigma);
    intercept_1 ~ normal(intercept_mu, intercept_sigma);
    calib_sigma ~ normal(0, 0.1);

    // Measurement priors
    od_per_biomass_mu ~ std_normal();
    od_per_biomass_sigma ~ normal(0, 0.1);

    // Likelihoods
    od_595nm_calib ~ normal(intercept_1[calib_idx] + slope_1[calib_idx] .* concentration, calib_sigma);
    od_per_biomass ~ normal(od_per_biomass_mu[meas_idx], od_per_biomass_sigma[meas_idx]);
}

generated quantities { 
    vector<lower=0>[N_meas] od_per_biomass_rep;
    vector<lower=0>[N_calib] od_595nm_calib_rep;
    //Posterior Predictive Checks
    for (i in 1:N_meas) {
        od_per_biomass_rep[i] = normal_rng(od_per_biomass_mu[meas_idx[i]], od_per_biomass_sigma[meas_idx[i]]);
    }
    for (i in 1:N_calib) { 
        od_595nm_calib_rep[i] = normal_rng(intercept_1[calib_idx[i]] + slope_1[calib_idx[i]] * concentration[i] , calib_sigma);
    } 
}